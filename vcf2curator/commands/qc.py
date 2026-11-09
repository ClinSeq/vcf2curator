"""``qc`` command: QC metric files for one case -> QC overview table.

Replaces the table output of the legacy ``QC_overview.R``. One run covers one
case: a tumor sample and an optional normal sample, each given its own Picard,
ContEst and (panel only) mSINGS files. The output has a header, the tumor row
and, if given, the normal row.

Behavioural notes vs. the legacy script (intentional, flagged for review):

* ``SAMP`` is the ``--tumor-name``/``--normal-name`` value as given, not
  shortened, and samples are never guessed from file or directory names.
* Rows are ordered tumor, then normal (legacy sorted by sample name).
* Picard files are read from their ``## METRICS CLASS`` section instead of
  skipping a fixed number of lines.
* Plots, the insert size histogram, flagstat input, directory scanning, and the
  small-design and rerun table variants are dropped.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Annotated, Any

import typer

from vcf2curator.commands._shared import QuietOpt, VerboseOpt, fail
from vcf2curator.curator.schema import QC_PANEL_COLUMNS, QC_WGS_COLUMNS
from vcf2curator.curator.writer import CuratorWriter, Row
from vcf2curator.exceptions import InputValidationError, Vcf2CuratorError
from vcf2curator.logging_config import configure_logging, get_logger
from vcf2curator.utils.picard import PicardMetrics, read_picard_metrics

logger = get_logger()

# Written for metrics that are absent or could not be computed (as R's write.table).
NA = "NA"

# mSINGS reports its summary as label/value pairs in the first lines.
_MSINGS_SUMMARY_LINES = 5


class QcMode(StrEnum):
    """Sequencing design, which decides the input files and output columns."""

    PANEL = "panel"
    WGS = "wgs"


@dataclass(frozen=True)
class SampleQcFiles:
    """The QC files for one sample of the case.

    Attributes:
        name: Sample name written to the SAMP column.
        metrics: Picard HsMetrics (panel) or WgsMetrics (WGS) file.
        insertsize: Picard InsertSize metrics file.
        contest: ContEst contamination file.
        markdup: Picard MarkDuplicates metrics file (panel only).
        msings: mSINGS ``MSI_Analysis.txt`` file (panel only, optional).
    """

    name: str
    metrics: Path
    insertsize: Path
    contest: Path
    markdup: Path | None = None
    msings: Path | None = None


# -- readers ----------------------------------------------------------------


def _value(metrics: PicardMetrics, column: str, path: Path) -> str:
    """Return a metric as written in the file, or NA if Picard left it empty."""
    if column not in metrics:
        raise InputValidationError(f"{path}: missing column '{column}'")
    value = metrics[column]
    return NA if value is None else value


def read_contamination(path: Path) -> str:
    """Return the ``contamination`` value from the first row of a ContEst file.

    Raises:
        InputValidationError: If the file cannot be read or lacks the column.
    """
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            first = next(csv.DictReader(handle, delimiter="\t"), None)
    except (OSError, UnicodeDecodeError) as exc:
        raise InputValidationError(f"Could not read ContEst file {path}: {exc}") from exc
    if first is None or "contamination" not in first:
        raise InputValidationError(
            f"{path}: expected a ContEst table with a 'contamination' column"
        )
    return first["contamination"].strip() or NA


def read_msing_score(path: Path) -> str:
    """Return ``msing_score`` from the label/value summary at the top of an mSINGS file.

    Raises:
        InputValidationError: If the file cannot be read or has no msing_score.
    """
    try:
        with path.open("r", encoding="utf-8") as handle:
            summary = [next(handle, "") for _ in range(_MSINGS_SUMMARY_LINES)]
    except (OSError, UnicodeDecodeError) as exc:
        raise InputValidationError(f"Could not read mSINGS file {path}: {exc}") from exc
    for line in summary:
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 2 and fields[0] == "msing_score":
            return fields[1].strip() or NA
    raise InputValidationError(f"{path}: no 'msing_score' line in the mSINGS summary")


def _ratio(numerator: str, denominator: str) -> str:
    """Return numerator/denominator to 15 significant digits, or NA."""
    try:
        value = float(numerator) / float(denominator)
    except (ValueError, ZeroDivisionError):
        return NA
    return f"{value:.15g}"


# -- rows -------------------------------------------------------------------


def build_qc_row(sample: SampleQcFiles, mode: QcMode) -> Row:
    """Build the QC overview row for one sample."""
    metrics = read_picard_metrics(sample.metrics)
    insert_size = read_picard_metrics(sample.insertsize)
    common = {
        "SAMP": sample.name,
        "FOLD_80_BASE_PENALTY": _value(metrics, "FOLD_80_BASE_PENALTY", sample.metrics),
        "contamination_%": read_contamination(sample.contest),
        "MEDIAN_INSERT_SIZE": _value(insert_size, "MEDIAN_INSERT_SIZE", sample.insertsize),
    }
    if mode is QcMode.WGS:
        return {**common, "MEAN_COVERAGE": _value(metrics, "MEAN_COVERAGE", sample.metrics)}

    assert sample.markdup is not None  # guaranteed by _sample_files
    markdup = read_picard_metrics(sample.markdup)
    return {
        **common,
        "MEAN_TARGET_COVERAGE": _value(metrics, "MEAN_TARGET_COVERAGE", sample.metrics),
        "FOLD_ENRICHMENT": _value(metrics, "FOLD_ENRICHMENT", sample.metrics),
        "dedupped_on_bait_rate": _ratio(
            _value(metrics, "ON_BAIT_BASES", sample.metrics),
            _value(metrics, "PF_BASES_ALIGNED", sample.metrics),
        ),
        "READ_PAIRS_EXAMINED": _value(markdup, "READ_PAIRS_EXAMINED", sample.markdup),
        "PERCENT_DUPLICATION": _value(markdup, "PERCENT_DUPLICATION", sample.markdup),
        "msing_score": read_msing_score(sample.msings) if sample.msings else NA,
    }


def columns_for_mode(mode: QcMode) -> tuple[str, ...]:
    """Return the output columns for a sequencing design."""
    return QC_WGS_COLUMNS if mode is QcMode.WGS else QC_PANEL_COLUMNS


# -- option handling --------------------------------------------------------


def _sample_files(
    role: str,
    name: str | None,
    hsmetrics: Path | None,
    wgsmetrics: Path | None,
    insertsize: Path | None,
    contest: Path | None,
    markdup: Path | None,
    msings: Path | None,
) -> tuple[SampleQcFiles, QcMode]:
    """Check one sample's options and bundle its files.

    Raises:
        InputValidationError: If a required option is missing or options conflict.
    """
    if (hsmetrics is None) == (wgsmetrics is None):
        raise InputValidationError(
            f"Give exactly one of --{role}-hsmetrics (panel) or --{role}-wgsmetrics (WGS)."
        )
    mode = QcMode.PANEL if hsmetrics is not None else QcMode.WGS
    required = {"name": name, "insertsize": insertsize, "contest": contest}
    if mode is QcMode.PANEL:
        required["markdup"] = markdup
    missing = [f"--{role}-{opt}" for opt, value in required.items() if value is None]
    if missing:
        raise InputValidationError(f"Missing required {mode.value} options: {', '.join(missing)}")
    if mode is QcMode.WGS and (markdup is not None or msings is not None):
        raise InputValidationError(
            f"--{role}-markdup and --{role}-msings are only used for panel data."
        )
    assert name is not None and insertsize is not None and contest is not None
    metrics = hsmetrics if hsmetrics is not None else wgsmetrics
    assert metrics is not None
    files = SampleQcFiles(name, metrics, insertsize, contest, markdup, msings)
    return files, mode


def _file_option(flag: str, help_text: str) -> Any:
    return typer.Option(flag, exists=True, dir_okay=False, help=help_text)


# Option types shared by the tumor and normal file options.
_FileOpt = Path | None


def qc(
    output: Annotated[Path, typer.Option("--output", "-o", help="Destination QC table.")],
    tumor_name: Annotated[
        str | None, typer.Option("--tumor-name", help="Tumor sample name (SAMP column).")
    ] = None,
    tumor_hsmetrics: Annotated[
        _FileOpt, _file_option("--tumor-hsmetrics", "Tumor Picard HsMetrics (panel).")
    ] = None,
    tumor_wgsmetrics: Annotated[
        _FileOpt, _file_option("--tumor-wgsmetrics", "Tumor Picard WgsMetrics (WGS).")
    ] = None,
    tumor_markdup: Annotated[
        _FileOpt, _file_option("--tumor-markdup", "Tumor Picard MarkDuplicates metrics (panel).")
    ] = None,
    tumor_insertsize: Annotated[
        _FileOpt, _file_option("--tumor-insertsize", "Tumor Picard InsertSize metrics.")
    ] = None,
    tumor_contest: Annotated[
        _FileOpt, _file_option("--tumor-contest", "Tumor ContEst contamination file.")
    ] = None,
    tumor_msings: Annotated[
        _FileOpt, _file_option("--tumor-msings", "Tumor mSINGS MSI_Analysis.txt (panel).")
    ] = None,
    normal_name: Annotated[
        str | None, typer.Option("--normal-name", help="Normal sample name (SAMP column).")
    ] = None,
    normal_hsmetrics: Annotated[
        _FileOpt, _file_option("--normal-hsmetrics", "Normal Picard HsMetrics (panel).")
    ] = None,
    normal_wgsmetrics: Annotated[
        _FileOpt, _file_option("--normal-wgsmetrics", "Normal Picard WgsMetrics (WGS).")
    ] = None,
    normal_markdup: Annotated[
        _FileOpt, _file_option("--normal-markdup", "Normal Picard MarkDuplicates metrics (panel).")
    ] = None,
    normal_insertsize: Annotated[
        _FileOpt, _file_option("--normal-insertsize", "Normal Picard InsertSize metrics.")
    ] = None,
    normal_contest: Annotated[
        _FileOpt, _file_option("--normal-contest", "Normal ContEst contamination file.")
    ] = None,
    normal_msings: Annotated[
        _FileOpt, _file_option("--normal-msings", "Normal mSINGS MSI_Analysis.txt (panel).")
    ] = None,
    dry_run: Annotated[
        bool, typer.Option("--dry-run", help="Parse and report counts without writing output.")
    ] = False,
    verbose: VerboseOpt = False,
    quiet: QuietOpt = False,
) -> None:
    """Combine one case's QC metric files (tumor, optional normal) into a QC table."""
    logger = configure_logging(verbose=verbose, quiet=quiet)
    normal_options = (
        normal_name,
        normal_hsmetrics,
        normal_wgsmetrics,
        normal_markdup,
        normal_insertsize,
        normal_contest,
        normal_msings,
    )

    try:
        tumor, mode = _sample_files(
            "tumor",
            tumor_name,
            tumor_hsmetrics,
            tumor_wgsmetrics,
            tumor_insertsize,
            tumor_contest,
            tumor_markdup,
            tumor_msings,
        )
        samples = [tumor]
        if any(opt is not None for opt in normal_options):
            normal, normal_mode = _sample_files(
                "normal",
                normal_name,
                normal_hsmetrics,
                normal_wgsmetrics,
                normal_insertsize,
                normal_contest,
                normal_markdup,
                normal_msings,
            )
            if normal_mode is not mode:
                raise InputValidationError(
                    f"Tumor is {mode.value} but normal is {normal_mode.value}; use the same design."
                )
            samples.append(normal)

        rows = [build_qc_row(sample, mode) for sample in samples]
        if dry_run:
            logger.info("Dry run: %d rows would be written to %s", len(rows), output)
            typer.echo(f"[dry-run] {len(rows)} rows")
            return
        count = CuratorWriter(columns_for_mode(mode)).write(rows, output)
        typer.echo(f"Wrote {count} rows to {output}")
    except Vcf2CuratorError as exc:
        fail(exc)
