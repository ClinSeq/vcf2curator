"""Reader for Picard metrics files (HsMetrics, WgsMetrics, MarkDuplicates, InsertSize).

A Picard metrics file has comment lines, then a ``## METRICS CLASS`` line, a
tab-separated header and one or more data rows, optionally followed by a
``## HISTOGRAM`` section. Only the first data row is used.
"""

from __future__ import annotations

from pathlib import Path

from vcf2curator.exceptions import InputValidationError

# Values Picard writes when a metric could not be computed.
_MISSING = frozenset({"", "?"})

# Column name -> value, with uncomputed metrics as None.
PicardMetrics = dict[str, str | None]


def read_picard_metrics(filepath: str | Path) -> PicardMetrics:
    """Return the first metrics row of a Picard metrics file.

    Args:
        filepath: Path to the Picard metrics file.

    Returns:
        Mapping of column name -> value; ``""`` and ``?`` become ``None``.

    Raises:
        InputValidationError: If the file cannot be read or has no metrics row.
    """
    path = Path(filepath)
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeDecodeError) as exc:
        raise InputValidationError(f"Could not read Picard metrics {path}: {exc}") from exc

    start = next((i for i, line in enumerate(lines) if line.startswith("## METRICS CLASS")), None)
    if start is None or start + 2 >= len(lines) or not lines[start + 2].strip():
        raise InputValidationError(f"{path}: no '## METRICS CLASS' section with a data row")

    header = lines[start + 1].split("\t")
    values = lines[start + 2].split("\t")
    values += [""] * (len(header) - len(values))
    return {
        col: (None if val.strip() in _MISSING else val.strip())
        for col, val in zip(header, values, strict=False)
    }
