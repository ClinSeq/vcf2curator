"""Tests for the shared Curator writer and its schema validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from vcf2curator.curator.schema import SOMATIC_COLUMNS
from vcf2curator.curator.writer import CuratorWriter
from vcf2curator.exceptions import SchemaError


def _complete_row() -> dict[str, object]:
    return {col: col.lower() for col in SOMATIC_COLUMNS}


def test_header_uses_schema_order() -> None:
    writer = CuratorWriter(SOMATIC_COLUMNS)
    assert writer.header().split("\t") == list(SOMATIC_COLUMNS)


def test_format_row_orders_and_serializes() -> None:
    writer = CuratorWriter(SOMATIC_COLUMNS)
    row = _complete_row()
    row["CHROM"] = 7
    row["N_VAF"] = None
    cells = writer.format_row(row).split("\t")
    assert cells[0] == "7"
    # None is rendered as the empty null token.
    assert cells[SOMATIC_COLUMNS.index("N_VAF")] == ""


def test_missing_column_raises_schema_error() -> None:
    writer = CuratorWriter(SOMATIC_COLUMNS)
    row = _complete_row()
    del row["GENE"]
    with pytest.raises(SchemaError, match="missing columns"):
        writer.format_row(row)


def test_extra_column_raises_schema_error() -> None:
    writer = CuratorWriter(SOMATIC_COLUMNS)
    row = _complete_row()
    row["SURPRISE"] = "x"
    with pytest.raises(SchemaError, match="unexpected columns"):
        writer.format_row(row)


def test_write_returns_row_count(tmp_path: Path) -> None:
    writer = CuratorWriter(SOMATIC_COLUMNS)
    out = tmp_path / "out.tsv"
    count = writer.write([_complete_row(), _complete_row()], out)
    assert count == 2
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "\t".join(SOMATIC_COLUMNS)
    assert len(lines) == 3  # header + 2 rows
