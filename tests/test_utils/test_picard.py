"""Tests for the Picard metrics reader."""

from __future__ import annotations

from pathlib import Path

import pytest

from vcf2curator.exceptions import InputValidationError
from vcf2curator.utils.picard import read_picard_metrics


def test_reads_first_row_only(fixtures_dir: Path) -> None:
    metrics = read_picard_metrics(fixtures_dir / "qc" / "tumor.markdup.txt")
    assert metrics["LIBRARY"] == "lib1"
    assert metrics["READ_PAIRS_EXAMINED"] == "45000000"


def test_question_mark_and_empty_become_none(fixtures_dir: Path) -> None:
    metrics = read_picard_metrics(fixtures_dir / "qc" / "normal.hsmetrics.txt")
    assert metrics["FOLD_80_BASE_PENALTY"] is None
    assert metrics["SAMPLE"] is None
    assert metrics["MEAN_TARGET_COVERAGE"] == "250.1"


def test_histogram_section_is_ignored(fixtures_dir: Path) -> None:
    metrics = read_picard_metrics(fixtures_dir / "qc" / "tumor.insertsize.txt")
    assert metrics["MEDIAN_INSERT_SIZE"] == "168"
    assert "insert_size" not in metrics


def test_short_row_is_padded(tmp_path: Path) -> None:
    path = tmp_path / "short.txt"
    path.write_text("## METRICS CLASS\tx\nA\tB\tC\n1\t2\n")
    assert read_picard_metrics(path) == {"A": "1", "B": "2", "C": None}


@pytest.mark.parametrize(
    "content",
    [
        "# no metrics section\n",
        "## METRICS CLASS\tx\nA\tB\n",  # header but no data row
        "## METRICS CLASS\tx\nA\tB\n\n",  # blank line instead of data
    ],
)
def test_missing_metrics_section(tmp_path: Path, content: str) -> None:
    path = tmp_path / "bad.txt"
    path.write_text(content)
    with pytest.raises(InputValidationError, match="METRICS CLASS"):
        read_picard_metrics(path)


def test_unreadable_file(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError):
        read_picard_metrics(tmp_path / "missing.txt")
