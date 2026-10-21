"""Shared, schema-validated writer for Curator output tables.

Every subcommand emits its table through this writer so that delimiter, column
order, header names, null handling and encoding are defined in exactly one place.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

from vcf2curator.exceptions import SchemaError
from vcf2curator.logging_config import get_logger

logger = get_logger()

Row = Mapping[str, Any]

DELIMITER = "\t"
ENCODING = "utf-8"
NULL_TOKEN = ""


class CuratorWriter:
    """Serialize rows to a Curator TSV against a fixed column schema.

    Args:
        columns: The ordered column names that make up the output header. Every
            row is emitted in exactly this order.
    """

    def __init__(self, columns: Sequence[str]) -> None:
        self.columns = tuple(columns)

    def _validate_row(self, row: Row) -> None:
        """Ensure a row contains exactly the schema columns.

        Raises:
            SchemaError: If the row is missing columns or has unexpected extras.
        """
        keys = set(row)
        expected = set(self.columns)
        missing = expected - keys
        extra = keys - expected
        if missing or extra:
            parts = []
            if missing:
                parts.append(f"missing columns {sorted(missing)}")
            if extra:
                parts.append(f"unexpected columns {sorted(extra)}")
            raise SchemaError("Row does not match Curator schema: " + "; ".join(parts))

    def _format_cell(self, value: Any) -> str:
        """Render a single cell, mapping null-ish values to the null token."""
        if value is None:
            return NULL_TOKEN
        return str(value)

    def format_row(self, row: Row) -> str:
        """Validate and render a single row as a delimited line (no newline)."""
        self._validate_row(row)
        return DELIMITER.join(self._format_cell(row[col]) for col in self.columns)

    def header(self) -> str:
        """Return the header line (no trailing newline)."""
        return DELIMITER.join(self.columns)

    def write(self, rows: Iterable[Row], output_path: str | Path) -> int:
        """Write the header and all rows to ``output_path``.

        Args:
            rows: An iterable of rows, each keyed by column name.
            output_path: Destination file path.

        Returns:
            The number of data rows written (excluding the header).
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        count = 0
        with output_path.open("w", encoding=ENCODING, newline="\n") as handle:
            handle.write(self.header() + "\n")
            for row in rows:
                handle.write(self.format_row(row) + "\n")
                count += 1
        logger.info("Wrote %d rows to %s", count, output_path)
        return count
