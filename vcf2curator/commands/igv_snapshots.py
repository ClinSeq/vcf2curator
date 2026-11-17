"""``igv-snapshots`` command: Curator table + IGV session -> PNG snapshots.

Replaces the legacy ``generate_igvsnapshots.py``. ``igv-snapshots snvs`` and
``igv-snapshots svs`` write an IGV batch script into ``--outdir`` and run IGV on
it, which saves one PNG per variant (one per breakpoint for SVs) there.

Behavioural notes vs. the legacy script (intentional, flagged for review):

* Columns are looked up by header name, not by position.
* SVs: rows with ``CURATOR == NO`` are skipped and every breakpoint uses the
  GRIDSS window (``END`` +/- 150 bp); the svcaller window is dropped, since the
  table no longer has a TOOL column.
* ``chr`` is added to the chromosome only when missing.
* Indels: ALT is split on ``,``; when the first ALT is shorter than REF
  (deletion) the view is centred at ``END + round(len(REF) / 2)``, otherwise at
  ``END``. Same result as the legacy code for real inputs.
* Rows giving an already-used snapshot name are skipped with a warning.
* Genome, IGV launcher and port are options; the preferences file is bundled
  and can be overridden.
* IGV is run without a shell; a failed run exits non-zero.
"""

from __future__ import annotations

import csv
import socket
import subprocess
from collections.abc import Callable, Iterable, Iterator, Sequence
from importlib.resources import as_file, files
from pathlib import Path
from typing import Annotated, Any

import typer

from vcf2curator.commands._shared import QuietOpt, VerboseOpt, fail
from vcf2curator.exceptions import IgvRunError, InputValidationError, Vcf2CuratorError
from vcf2curator.logging_config import configure_logging, get_logger

logger = get_logger()

SNV_BATCH_NAME = "igv_batch_snv.sh"
SV_BATCH_NAME = "igv_batch_sv.sh"
SV_FLANK = 150

SNV_REQUIRED_COLUMNS: tuple[str, ...] = ("CHROM", "START", "END", "REF", "ALT")
SV_REQUIRED_COLUMNS: tuple[str, ...] = (
    "CHROM_A",
    "START_A",
    "END_A",
    "CHROM_B",
    "START_B",
    "END_B",
    "CURATOR",
)

_PREFS_RESOURCE = files("vcf2curator.resources.igv").joinpath("igvsnapshot.prefs")

Row = dict[str, str]

igv_snapshots_app = typer.Typer(
    no_args_is_help=True,
    help="Take IGV snapshots of the variants in a Curator table.",
)


# -- table ------------------------------------------------------------------


def read_table(path: Path, required: Sequence[str]) -> list[Row]:
    """Read a tab-separated Curator table into rows keyed by column name.

    Args:
        path: Table file.
        required: Columns that must be present in the header.

    Returns:
        One dict per data row.

    Raises:
        InputValidationError: If the file is unreadable, has no header or lacks
            a required column.
    """
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = reader.fieldnames
            if not header:
                raise InputValidationError(f"{path}: no header line")
            missing = [col for col in required if col not in header]
            if missing:
                raise InputValidationError(f"{path}: missing columns {missing}")
            return [{k: (v or "") for k, v in row.items() if k is not None} for row in reader]
    except (OSError, UnicodeDecodeError) as exc:
        raise InputValidationError(f"Could not read {path}: {exc}") from exc


def _int(row: Row, column: str, line: int) -> int:
    try:
        return int(row[column])
    except ValueError as exc:
        raise InputValidationError(
            f"Row {line}: {column} must be an integer, got {row[column]!r}"
        ) from exc


# -- commands ---------------------------------------------------------------


def ucsc_chrom(chrom: str) -> str:
    """Return ``chrom`` with a ``chr`` prefix, adding it only when missing."""
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def snapshot_name(chrom: str, start: int, end: int) -> str:
    """Return the legacy snapshot file name for a locus."""
    return f"{ucsc_chrom(chrom)}_{start}-{end}.png"


def snv_goto_pos(ref: str, alt: str, end: int) -> int:
    """Return the position IGV is centred on for an SNV/indel.

    Deletions (first ALT shorter than REF) are centred inside the deleted
    bases; everything else at ``end``.
    """
    first_alt = alt.split(",")[0]
    if len(first_alt) < len(ref):
        return end + round(len(ref) / 2)
    return end


def _snapshot(locus: str, name: str) -> list[str]:
    return [f"goto {locus}", "sort base", f"snapshot {name}"]


class _Dedup:
    """Tracks snapshot names and warns about repeats."""

    def __init__(self) -> None:
        self._seen: set[str] = set()

    def is_new(self, name: str) -> bool:
        if name in self._seen:
            logger.warning("Duplicate snapshot %s; skipping", name)
            return False
        self._seen.add(name)
        return True


def snv_commands(rows: Iterable[Row]) -> list[str]:
    """Build the batch commands for an SNV/indel table.

    Raises:
        InputValidationError: If START or END is not an integer.
    """
    cmds: list[str] = []
    dedup = _Dedup()
    for line, row in enumerate(rows, start=2):
        start, end = _int(row, "START", line), _int(row, "END", line)
        name = snapshot_name(row["CHROM"], start, end)
        if not dedup.is_new(name):
            continue
        pos = snv_goto_pos(row["REF"], row["ALT"], end)
        cmds.extend(_snapshot(f"{ucsc_chrom(row['CHROM'])}:{pos}", name))
    return cmds


def _sv_breakpoints(row: Row, line: int) -> Iterator[tuple[str, int, int]]:
    for side in ("A", "B"):
        yield (
            row[f"CHROM_{side}"],
            _int(row, f"START_{side}", line),
            _int(row, f"END_{side}", line),
        )


def sv_commands(rows: Iterable[Row]) -> list[str]:
    """Build the batch commands for an SV table (CURATOR=NO rows skipped).

    Raises:
        InputValidationError: If a START/END column is not an integer.
    """
    cmds: list[str] = ["squish tb_simpleRepeat"]
    dedup = _Dedup()
    for line, row in enumerate(rows, start=2):
        if row["CURATOR"] == "NO":
            continue
        for chrom, start, end in _sv_breakpoints(row, line):
            name = snapshot_name(chrom, start, end)
            if not dedup.is_new(name):
                continue
            locus = f"{ucsc_chrom(chrom)}:{end - SV_FLANK}-{end + SV_FLANK}"
            cmds.extend(_snapshot(locus, name))
    return cmds


def count_snapshots(cmds: Sequence[str]) -> int:
    """Return the number of ``snapshot`` commands in ``cmds``."""
    return sum(1 for cmd in cmds if cmd.startswith("snapshot "))


def batch_script(body: Sequence[str], *, genome: str, outdir: Path, session: Path) -> str:
    """Wrap per-variant commands with the session setup and ``exit``."""
    lines = [
        "new",
        f"genome {genome}",
        f"snapshotDirectory {outdir.resolve()}",
        f"load {session.resolve()}",
        *body,
        "exit",
    ]
    return "\n".join(lines) + "\n"


# -- IGV --------------------------------------------------------------------


def free_port() -> int:
    """Return a TCP port that is free right now."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(("127.0.0.1", 0))
        port: int = sock.getsockname()[1]
        return port


def igv_command(igv: str, batch: Path, prefs: Path, port: int, outdir: Path) -> list[str]:
    """Return the IGV command line."""
    return [
        igv,
        "-b",
        str(batch),
        "--preferences",
        str(prefs),
        "-p",
        str(port),
        "--igvDirectory",
        str(outdir),
    ]


def run_igv(cmd: Sequence[str]) -> None:
    """Run IGV.

    Raises:
        IgvRunError: If IGV cannot be started or exits non-zero.
    """
    logger.info("Running %s", " ".join(cmd))
    try:
        result = subprocess.run(list(cmd), check=False)
    except OSError as exc:
        raise IgvRunError(f"Could not start IGV ({cmd[0]}): {exc}") from exc
    if result.returncode != 0:
        raise IgvRunError(f"IGV exited with code {result.returncode}")


# -- CLI --------------------------------------------------------------------

InputOpt = Annotated[
    Path,
    typer.Option("--input", "-i", exists=True, dir_okay=False, help="Curator table (TSV)."),
]
SessionOpt = Annotated[
    Path,
    typer.Option("--session", exists=True, dir_okay=False, help="IGV session XML."),
]
OutdirOpt = Annotated[
    Path,
    typer.Option("--outdir", "-o", file_okay=False, help="Directory for the snapshots."),
]
PrefsOpt = Annotated[
    Path | None,
    typer.Option(
        "--preferences",
        exists=True,
        dir_okay=False,
        help="IGV preferences file (default: bundled igvsnapshot.prefs).",
    ),
]
GenomeOpt = Annotated[str, typer.Option("--genome", help="IGV genome ID.")]
IgvOpt = Annotated[str, typer.Option("--igv", help="IGV launcher to run.")]
PortOpt = Annotated[
    int | None,
    typer.Option("--port", min=1, max=65535, help="IGV batch port (default: a free port)."),
]
DryRunOpt = Annotated[
    bool, typer.Option("--dry-run", help="Log the batch script without writing it or running IGV.")
]


def _run(
    *,
    body: list[str],
    batch_name: str,
    session: Path,
    outdir: Path,
    preferences: Path | None,
    genome: str,
    igv: str,
    port: int | None,
    dry_run: bool,
) -> None:
    """Write the batch script and run IGV on it."""
    script = batch_script(body, genome=genome, outdir=outdir, session=session)
    batch = outdir.resolve() / batch_name
    n_snapshots = count_snapshots(body)
    port = port if port is not None else free_port()

    with as_file(_PREFS_RESOURCE) as bundled_prefs:
        prefs = (preferences or bundled_prefs).resolve()
        cmd = igv_command(igv, batch, prefs, port, outdir.resolve())

        if dry_run:
            logger.info(
                "Dry run: would write %s with %d snapshots:\n%s", batch, n_snapshots, script
            )
            logger.info("Dry run: would run %s", " ".join(cmd))
            return

        outdir.mkdir(parents=True, exist_ok=True)
        batch.write_text(script, encoding="utf-8")
        logger.info("Wrote %s (%d snapshots)", batch, n_snapshots)
        if n_snapshots == 0:
            logger.warning("No variants to snapshot; not running IGV")
            return
        run_igv(cmd)
    logger.info("Snapshots written to %s", outdir)


def _execute(
    build: Callable[[list[Row]], list[str]],
    input_path: Path,
    required: Sequence[str],
    **kwargs: Any,
) -> None:
    """Read the table, build the commands and run; domain errors exit non-zero."""
    try:
        body = build(read_table(input_path, required))
        _run(body=body, **kwargs)
    except Vcf2CuratorError as exc:
        fail(exc)


@igv_snapshots_app.command("snvs")
def snvs(
    input_path: InputOpt,
    session: SessionOpt,
    outdir: OutdirOpt,
    preferences: PrefsOpt = None,
    genome: GenomeOpt = "hg19",
    igv: IgvOpt = "igv.sh",
    port: PortOpt = None,
    dry_run: DryRunOpt = False,
    verbose: VerboseOpt = False,
    quiet: QuietOpt = False,
) -> None:
    """Snapshot every variant in an SNV/indel Curator table."""
    configure_logging(verbose=verbose, quiet=quiet)
    _execute(
        snv_commands,
        input_path,
        SNV_REQUIRED_COLUMNS,
        batch_name=SNV_BATCH_NAME,
        session=session,
        outdir=outdir,
        preferences=preferences,
        genome=genome,
        igv=igv,
        port=port,
        dry_run=dry_run,
    )


@igv_snapshots_app.command("svs")
def svs(
    input_path: InputOpt,
    session: SessionOpt,
    outdir: OutdirOpt,
    preferences: PrefsOpt = None,
    genome: GenomeOpt = "hg19",
    igv: IgvOpt = "igv.sh",
    port: PortOpt = None,
    dry_run: DryRunOpt = False,
    verbose: VerboseOpt = False,
    quiet: QuietOpt = False,
) -> None:
    """Snapshot both breakpoints of every CURATOR=YES SV in a Curator table."""
    configure_logging(verbose=verbose, quiet=quiet)
    _execute(
        sv_commands,
        input_path,
        SV_REQUIRED_COLUMNS,
        batch_name=SV_BATCH_NAME,
        session=session,
        outdir=outdir,
        preferences=preferences,
        genome=genome,
        igv=igv,
        port=port,
        dry_run=dry_run,
    )
