"""Tests for BED loading and interval overlap helpers."""

from __future__ import annotations

from pathlib import Path

import pytest

from vcf2curator.annotations.intervals import load_bed, near_any, overlapping_names, overlaps_any
from vcf2curator.exceptions import MissingAnnotationError


def test_load_bed_groups_by_chrom(fixtures_dir: Path) -> None:
    genes = load_bed(fixtures_dir / "sv_genes.bed", require_name=True)
    assert sorted(genes) == ["1", "7"]
    assert [i.name for i in genes["1"]] == ["BRCA1", "GENEX"]


def test_load_bed_requires_name_column(fixtures_dir: Path) -> None:
    with pytest.raises(MissingAnnotationError, match="at least 4"):
        load_bed(fixtures_dir / "sv_targets.bed", require_name=True)


def test_load_bed_rejects_non_integer_coords(tmp_path: Path) -> None:
    bad = tmp_path / "bad.bed"
    bad.write_text("1\tx\t10\n")
    with pytest.raises(MissingAnnotationError, match="integers"):
        load_bed(bad)


def test_load_bed_missing_file(tmp_path: Path) -> None:
    with pytest.raises(MissingAnnotationError):
        load_bed(tmp_path / "missing.bed")


def test_overlaps_are_half_open(fixtures_dir: Path) -> None:
    genes = load_bed(fixtures_dir / "sv_genes.bed", require_name=True)
    assert overlapping_names(genes, "1", 999, 1000) == ["BRCA1"]
    assert overlapping_names(genes, "1", 1100, 1200) == []
    assert overlapping_names(genes, "1", 1000, 2000) == ["BRCA1", "GENEX"]
    assert overlaps_any(genes, "7", 599, 600)
    assert not overlaps_any(genes, "7", 600, 700)
    assert not overlaps_any(genes, "2", 0, 10_000)


def test_near_any_uses_padding(fixtures_dir: Path) -> None:
    targets = load_bed(fixtures_dir / "sv_targets.bed")
    assert near_any(targets, "1", (1200,), padding=150)
    assert not near_any(targets, "1", (1201,), padding=150)
    assert near_any(targets, "1", (3000, 800), padding=150)
