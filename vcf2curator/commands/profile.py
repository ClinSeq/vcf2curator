"""``profile`` command: one case's vcf2curator outputs -> a Curator YAML profile.

The profile is the Curator app's single entry point for a case: metadata plus
the paths of every output (SNV/SV/QC tables, IGV sessions and snapshots,
Frankenplot and MultiQC reports). Files are referenced, not embedded.

Layout (``schema_version: 1``; sections and keys with no input are omitted)::

    schema_version: 1
    case_id: ...
    project: ...
    cancer_type: ...
    samples: {tumor: ..., normal: ...}
    variants:
      snvs: {somatic: ..., germline: ...}
      svs: {somatic: ..., germline: ..., evidence_bam: ...}
    qc: {table: ..., multiqc: ...}
    igv:
      sessions: {snvs: ..., svs: ...}
      snapshots: {snvs: ..., svs: ...}
    cnv: {frankenplot_html: ..., frankenplot_images: [...]}
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from dataclasses import dataclass, field, fields
from enum import StrEnum
from pathlib import Path
from typing import Annotated, Any

import typer
import yaml

from vcf2curator.commands._shared import QuietOpt, VerboseOpt, fail
from vcf2curator.curator.schema import (
    GERMLINE_COLUMNS,
    QC_PANEL_COLUMNS,
    QC_WGS_COLUMNS,
    SOMATIC_COLUMNS,
    SV_COLUMNS,
)
from vcf2curator.curator.writer import DELIMITER, ENCODING
from vcf2curator.exceptions import InputValidationError, Vcf2CuratorError
from vcf2curator.logging_config import configure_logging, get_logger

logger = get_logger()

SCHEMA_VERSION = 1


class PathStyle(StrEnum):
    """How file paths are written into the profile."""

    ABSOLUTE = "absolute"
    RELATIVE = "relative"


@dataclass(frozen=True)
class ProfileInputs:
    """Metadata and output files of one case.

    Attributes:
        case_id: Case identifier.
        tumor_name: Tumor sample name (as in the QC table's SAMP column).
        snvs_somatic: Somatic SNV/indel Curator table.
        qc: QC overview table.
        project: Project name.
        cancer_type: Cancer type (free text).
        normal_name: Normal sample name; ``None`` for tumor-only cases.
        snvs_germline: Germline SNV/indel Curator table.
        svs_somatic: Somatic SV Curator table.
        svs_germline: Germline SV Curator table.
        sv_evidence_bam: GRIDSS evidence BAM.
        multiqc: MultiQC HTML report.
        igv_session_snvs: IGV session XML for SNVs.
        igv_session_svs: IGV session XML for SVs.
        igv_snapshots_snvs: Directory of SNV IGV snapshots.
        igv_snapshots_svs: Directory of SV IGV snapshots.
        frankenplot_html: Frankenplot HTML report.
        frankenplot_images: Frankenplot images.
    """

    case_id: str
    tumor_name: str
    snvs_somatic: Path
    qc: Path
    project: str | None = None
    cancer_type: str | None = None
    normal_name: str | None = None
    snvs_germline: Path | None = None
    svs_somatic: Path | None = None
    svs_germline: Path | None = None
    sv_evidence_bam: Path | None = None
    multiqc: Path | None = None
    igv_session_snvs: Path | None = None
    igv_session_svs: Path | None = None
    igv_snapshots_snvs: Path | None = None
    igv_snapshots_svs: Path | None = None
    frankenplot_html: Path | None = None
    frankenplot_images: tuple[Path, ...] = field(default_factory=tuple)


# Tables whose header must match a Curator schema, by ProfileInputs attribute.
_TABLE_SCHEMAS: dict[str, tuple[tuple[str, ...], ...]] = {
    "snvs_somatic": (SOMATIC_COLUMNS,),
    "snvs_germline": (GERMLINE_COLUMNS,),
    "svs_somatic": (SV_COLUMNS,),
    "svs_germline": (SV_COLUMNS,),
    "qc": (QC_PANEL_COLUMNS, QC_WGS_COLUMNS),
}

_SNAPSHOT_DIRS = ("igv_snapshots_snvs", "igv_snapshots_svs")


# -- validation -------------------------------------------------------------


def _read_lines(path: Path) -> list[list[str]]:
    try:
        text = path.read_text(encoding=ENCODING)
    except (OSError, UnicodeDecodeError) as exc:
        raise InputValidationError(f"Could not read {path}: {exc}") from exc
    return [line.split(DELIMITER) for line in text.splitlines() if line]


def check_table_header(path: Path, allowed: Sequence[tuple[str, ...]]) -> tuple[str, ...]:
    """Check that a table's header equals one of the allowed Curator schemas.

    Args:
        path: Tab-separated table.
        allowed: Accepted column tuples.

    Returns:
        The matching column tuple.

    Raises:
        InputValidationError: If the file is empty or its header matches none.
    """
    lines = _read_lines(path)
    if not lines:
        raise InputValidationError(f"{path}: empty table")
    header = tuple(lines[0])
    for columns in allowed:
        if header == columns:
            return columns
    expected = " or ".join(
        f"[{', '.join(cols[:3])}, ...] ({len(cols)} columns)" for cols in allowed
    )
    raise InputValidationError(
        f"{path}: header does not match the expected Curator table {expected}; "
        f"got [{', '.join(header[:3])}, ...] ({len(header)} columns)"
    )


def check_qc_samples(path: Path, names: Sequence[str]) -> None:
    """Check that every sample name appears in the QC table's SAMP column.

    Raises:
        InputValidationError: If a name is missing.
    """
    lines = _read_lines(path)
    samp = lines[0].index("SAMP")
    present = {row[samp] for row in lines[1:] if len(row) > samp}
    missing = [name for name in names if name not in present]
    if missing:
        raise InputValidationError(
            f"{path}: sample(s) {missing} not in SAMP column (found {sorted(present)})"
        )


def _check_not_empty(inputs: ProfileInputs) -> None:
    for f in fields(inputs):
        value = getattr(inputs, f.name)
        paths = value if isinstance(value, tuple) else (value,)
        for path in paths:
            if isinstance(path, Path) and path.is_file() and path.stat().st_size == 0:
                raise InputValidationError(f"Empty file for --{f.name.replace('_', '-')}: {path}")


def _warn_empty_snapshot_dirs(inputs: ProfileInputs) -> None:
    for name in _SNAPSHOT_DIRS:
        directory: Path | None = getattr(inputs, name)
        if directory is not None and not any(directory.glob("*.png")):
            logger.warning("No PNG snapshots in %s", directory)


def _warn_missing_bam_index(bam: Path | None) -> None:
    if bam is None:
        return
    if not (bam.with_name(bam.name + ".bai").exists() or bam.with_suffix(".bai").exists()):
        logger.warning("No .bai index next to %s", bam)


def validate_inputs(inputs: ProfileInputs) -> None:
    """Validate a case's files before writing the profile.

    Raises:
        InputValidationError: On an empty file, a wrong table header or QC
            sample names that don't match.
    """
    _check_not_empty(inputs)
    for name, allowed in _TABLE_SCHEMAS.items():
        path: Path | None = getattr(inputs, name)
        if path is not None:
            check_table_header(path, allowed)
    names = [inputs.tumor_name] + ([inputs.normal_name] if inputs.normal_name else [])
    check_qc_samples(inputs.qc, names)
    _warn_empty_snapshot_dirs(inputs)
    _warn_missing_bam_index(inputs.sv_evidence_bam)


# -- rendering --------------------------------------------------------------


def format_path(path: Path, base: Path, style: PathStyle) -> str:
    """Render a path as absolute or relative to ``base``.

    Symlinks are not followed: Nextflow publishes outputs as symlinks whose
    targets live in work directories that may be cleaned up.
    """
    absolute = os.path.abspath(path)
    if style is PathStyle.RELATIVE:
        return os.path.relpath(absolute, os.path.abspath(base))
    return absolute


def _prune(mapping: dict[str, Any]) -> dict[str, Any]:
    """Drop ``None`` values and empty sections, recursively."""
    pruned: dict[str, Any] = {}
    for key, value in mapping.items():
        if isinstance(value, dict):
            value = _prune(value)
        if value is None or value == {} or value == []:
            continue
        pruned[key] = value
    return pruned


def build_profile(inputs: ProfileInputs, *, base: Path, style: PathStyle) -> dict[str, Any]:
    """Build the profile mapping in the documented key order.

    Args:
        inputs: Case metadata and files.
        base: Directory the profile is written to (base for relative paths).
        style: Absolute or relative paths.
    """

    def p(path: Path | None) -> str | None:
        return None if path is None else format_path(path, base, style)

    return _prune(
        {
            "schema_version": SCHEMA_VERSION,
            "case_id": inputs.case_id,
            "project": inputs.project,
            "cancer_type": inputs.cancer_type,
            "samples": {"tumor": inputs.tumor_name, "normal": inputs.normal_name},
            "variants": {
                "snvs": {"somatic": p(inputs.snvs_somatic), "germline": p(inputs.snvs_germline)},
                "svs": {
                    "somatic": p(inputs.svs_somatic),
                    "germline": p(inputs.svs_germline),
                    "evidence_bam": p(inputs.sv_evidence_bam),
                },
            },
            "qc": {"table": p(inputs.qc), "multiqc": p(inputs.multiqc)},
            "igv": {
                "sessions": {"snvs": p(inputs.igv_session_snvs), "svs": p(inputs.igv_session_svs)},
                "snapshots": {
                    "snvs": p(inputs.igv_snapshots_snvs),
                    "svs": p(inputs.igv_snapshots_svs),
                },
            },
            "cnv": {
                "frankenplot_html": p(inputs.frankenplot_html),
                "frankenplot_images": [p(img) for img in inputs.frankenplot_images],
            },
        }
    )


def render_yaml(profile: dict[str, Any]) -> str:
    """Serialize a profile to YAML, keeping key order."""
    return yaml.safe_dump(profile, sort_keys=False, default_flow_style=False)


# -- CLI --------------------------------------------------------------------


def _file(flag: str, help_text: str) -> Any:
    return typer.Option(flag, exists=True, dir_okay=False, help=help_text)


def _dir(flag: str, help_text: str) -> Any:
    return typer.Option(flag, exists=True, file_okay=False, help=help_text)


def profile(
    case_id: Annotated[str, typer.Option("--case-id", help="Case identifier.")],
    tumor_name: Annotated[str, typer.Option("--tumor-name", help="Tumor sample name.")],
    snvs_somatic: Annotated[Path, _file("--snvs-somatic", "Somatic SNV Curator table.")],
    qc: Annotated[Path, _file("--qc", "QC overview table.")],
    output: Annotated[
        Path, typer.Option("--output", "-o", dir_okay=False, help="Destination YAML profile.")
    ],
    project: Annotated[str | None, typer.Option("--project", help="Project name.")] = None,
    cancer_type: Annotated[
        str | None, typer.Option("--cancer-type", help="Cancer type (free text).")
    ] = None,
    normal_name: Annotated[
        str | None, typer.Option("--normal-name", help="Normal sample name (omit if tumor-only).")
    ] = None,
    snvs_germline: Annotated[
        Path | None, _file("--snvs-germline", "Germline SNV Curator table.")
    ] = None,
    svs_somatic: Annotated[Path | None, _file("--svs-somatic", "Somatic SV Curator table.")] = None,
    svs_germline: Annotated[
        Path | None, _file("--svs-germline", "Germline SV Curator table.")
    ] = None,
    sv_evidence_bam: Annotated[
        Path | None, _file("--sv-evidence-bam", "GRIDSS evidence BAM.")
    ] = None,
    multiqc: Annotated[Path | None, _file("--multiqc", "MultiQC HTML report.")] = None,
    igv_session_snvs: Annotated[
        Path | None, _file("--igv-session-snvs", "IGV session XML for SNVs.")
    ] = None,
    igv_session_svs: Annotated[
        Path | None, _file("--igv-session-svs", "IGV session XML for SVs.")
    ] = None,
    igv_snapshots_snvs: Annotated[
        Path | None, _dir("--igv-snapshots-snvs", "Directory of SNV IGV snapshots.")
    ] = None,
    igv_snapshots_svs: Annotated[
        Path | None, _dir("--igv-snapshots-svs", "Directory of SV IGV snapshots.")
    ] = None,
    frankenplot_html: Annotated[
        Path | None, _file("--frankenplot-html", "Frankenplot HTML report.")
    ] = None,
    frankenplot_image: Annotated[
        list[Path] | None,
        _file("--frankenplot-image", "Frankenplot image; repeat for several."),
    ] = None,
    path_style: Annotated[
        PathStyle,
        typer.Option(
            "--path-style", help="Write paths absolute or relative to the profile's directory."
        ),
    ] = PathStyle.ABSOLUTE,
    dry_run: Annotated[
        bool, typer.Option("--dry-run", help="Validate and log the profile without writing it.")
    ] = False,
    verbose: VerboseOpt = False,
    quiet: QuietOpt = False,
) -> None:
    """Write a Curator YAML profile listing one case's outputs."""
    configure_logging(verbose=verbose, quiet=quiet)
    inputs = ProfileInputs(
        case_id=case_id,
        tumor_name=tumor_name,
        snvs_somatic=snvs_somatic,
        qc=qc,
        project=project,
        cancer_type=cancer_type,
        normal_name=normal_name,
        snvs_germline=snvs_germline,
        svs_somatic=svs_somatic,
        svs_germline=svs_germline,
        sv_evidence_bam=sv_evidence_bam,
        multiqc=multiqc,
        igv_session_snvs=igv_session_snvs,
        igv_session_svs=igv_session_svs,
        igv_snapshots_snvs=igv_snapshots_snvs,
        igv_snapshots_svs=igv_snapshots_svs,
        frankenplot_html=frankenplot_html,
        frankenplot_images=tuple(frankenplot_image or ()),
    )
    try:
        validate_inputs(inputs)
    except Vcf2CuratorError as exc:
        fail(exc)

    text = render_yaml(build_profile(inputs, base=output.parent, style=path_style))
    if dry_run:
        logger.info("Dry run: would write %s:\n%s", output, text)
        return

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text, encoding=ENCODING)
    logger.info("Wrote %s", output)
