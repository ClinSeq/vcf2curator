"""Tests for the QC overview table."""

from __future__ import annotations

from pathlib import Path

import pytest

from vcf2curator.commands.qc import (
    QcMode,
    SampleQcFiles,
    build_qc_row,
    columns_for_mode,
    read_contamination,
    read_msing_score,
)
from vcf2curator.curator.schema import QC_PANEL_COLUMNS, QC_WGS_COLUMNS
from vcf2curator.exceptions import InputValidationError


@pytest.fixture
def qc_dir(fixtures_dir: Path) -> Path:
    return fixtures_dir / "qc"


def _panel_tumor(qc_dir: Path, msings: bool = True) -> SampleQcFiles:
    return SampleQcFiles(
        name="P-0001-T-01",
        metrics=qc_dir / "tumor.hsmetrics.txt",
        insertsize=qc_dir / "tumor.insertsize.txt",
        contest=qc_dir / "tumor.contest.txt",
        markdup=qc_dir / "tumor.markdup.txt",
        msings=qc_dir / "tumor.MSI_Analysis.txt" if msings else None,
    )


def test_panel_row(qc_dir: Path) -> None:
    row = build_qc_row(_panel_tumor(qc_dir), QcMode.PANEL)
    assert set(row) == set(QC_PANEL_COLUMNS)
    assert dict(row) == {
        "SAMP": "P-0001-T-01",
        "MEAN_TARGET_COVERAGE": "512.34",
        "FOLD_ENRICHMENT": "1210.5",
        "dedupped_on_bait_rate": "0.75",
        "FOLD_80_BASE_PENALTY": "1.45",
        "READ_PAIRS_EXAMINED": "45000000",
        "PERCENT_DUPLICATION": "0.4123",
        "contamination_%": "0.3",
        "MEDIAN_INSERT_SIZE": "168",
        "msing_score": "0.0500",
    }


def test_panel_row_without_msings_and_uncomputed_metric(qc_dir: Path) -> None:
    normal = SampleQcFiles(
        name="P-0001-N-01",
        metrics=qc_dir / "normal.hsmetrics.txt",
        insertsize=qc_dir / "normal.insertsize.txt",
        contest=qc_dir / "normal.contest.txt",
        markdup=qc_dir / "normal.markdup.txt",
    )
    row = build_qc_row(normal, QcMode.PANEL)
    assert row["msing_score"] == "NA"
    assert row["FOLD_80_BASE_PENALTY"] == "NA"  # "?" in the HsMetrics file
    assert row["dedupped_on_bait_rate"] == "0.6"


def test_wgs_row(qc_dir: Path) -> None:
    sample = SampleQcFiles(
        name="P-0001-T-01",
        metrics=qc_dir / "tumor.wgsmetrics.txt",
        insertsize=qc_dir / "tumor.insertsize.txt",
        contest=qc_dir / "tumor.contest.txt",
    )
    row = build_qc_row(sample, QcMode.WGS)
    assert set(row) == set(QC_WGS_COLUMNS)
    assert (row["MEAN_COVERAGE"], row["FOLD_80_BASE_PENALTY"]) == ("61.234", "1.21")


def test_columns_for_mode() -> None:
    assert columns_for_mode(QcMode.PANEL) == QC_PANEL_COLUMNS
    assert columns_for_mode(QcMode.WGS) == QC_WGS_COLUMNS


def test_hsmetrics_given_as_wgs_is_rejected(qc_dir: Path) -> None:
    sample = SampleQcFiles(
        name="x",
        metrics=qc_dir / "tumor.hsmetrics.txt",
        insertsize=qc_dir / "tumor.insertsize.txt",
        contest=qc_dir / "tumor.contest.txt",
    )
    with pytest.raises(InputValidationError, match="MEAN_COVERAGE"):
        build_qc_row(sample, QcMode.WGS)


def test_on_bait_ratio_na_when_denominator_zero(qc_dir: Path, tmp_path: Path) -> None:
    text = (qc_dir / "tumor.hsmetrics.txt").read_text().replace("800000000", "0")
    hs = tmp_path / "zero.hsmetrics.txt"
    hs.write_text(text)
    sample = SampleQcFiles(
        name="x",
        metrics=hs,
        insertsize=qc_dir / "tumor.insertsize.txt",
        contest=qc_dir / "tumor.contest.txt",
        markdup=qc_dir / "tumor.markdup.txt",
    )
    assert build_qc_row(sample, QcMode.PANEL)["dedupped_on_bait_rate"] == "NA"


def test_contest_without_contamination_column(tmp_path: Path) -> None:
    path = tmp_path / "bad.contest.txt"
    path.write_text("name\tvalue\nMETA\t0.1\n")
    with pytest.raises(InputValidationError, match="contamination"):
        read_contamination(path)


def test_msings_without_score(tmp_path: Path) -> None:
    path = tmp_path / "bad.MSI_Analysis.txt"
    path.write_text("Position\tx\nunstable_loci\t1\n")
    with pytest.raises(InputValidationError, match="msing_score"):
        read_msing_score(path)
