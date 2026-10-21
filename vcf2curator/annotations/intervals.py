"""BED interval loading and overlap checks for SV gene/exon/target annotation."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from vcf2curator.exceptions import MissingAnnotationError


@dataclass(frozen=True)
class Interval:
    """A 0-based, half-open BED interval with an optional name."""

    start: int
    end: int
    name: str = ""


# Intervals grouped by chromosome, in file order.
IntervalIndex = dict[str, list[Interval]]


def load_bed(filepath: str | Path, *, require_name: bool = False) -> IntervalIndex:
    """Load a BED file into intervals grouped by chromosome.

    Blank lines and ``#``/``track``/``browser`` header lines are skipped.

    Args:
        filepath: Path to the BED file.
        require_name: Fail if a line has no 4th (name) column.

    Returns:
        Mapping of chromosome -> intervals.

    Raises:
        MissingAnnotationError: If the file cannot be read or a line is malformed.
    """
    path = Path(filepath)
    index: IntervalIndex = defaultdict(list)
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line_no, line in enumerate(handle, start=1):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3 or (require_name and len(fields) < 4):
                    expected = 4 if require_name else 3
                    raise MissingAnnotationError(
                        f"{path}:{line_no}: expected at least {expected} tab-separated columns"
                    )
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError as exc:
                    raise MissingAnnotationError(
                        f"{path}:{line_no}: start/end must be integers"
                    ) from exc
                name = fields[3] if len(fields) > 3 else ""
                index[fields[0]].append(Interval(start, end, name))
    except OSError as exc:
        raise MissingAnnotationError(f"Could not read BED file {path}: {exc}") from exc
    return dict(index)


def overlapping_names(index: IntervalIndex, chrom: str, start: int, end: int) -> list[str]:
    """Return the unique names of intervals overlapping ``[start, end)``, in file order."""
    names: list[str] = []
    for interval in index.get(chrom, []):
        if interval.start < end and interval.end > start and interval.name not in names:
            names.append(interval.name)
    return names


def overlaps_any(index: IntervalIndex, chrom: str, start: int, end: int) -> bool:
    """Return True if any interval overlaps ``[start, end)``."""
    return any(i.start < end and i.end > start for i in index.get(chrom, []))


def near_any(index: IntervalIndex, chrom: str, positions: tuple[int, ...], padding: int) -> bool:
    """Return True if any position lies within ``padding`` bp of an interval."""
    return any(
        i.start - padding <= pos <= i.end + padding
        for i in index.get(chrom, [])
        for pos in positions
    )
