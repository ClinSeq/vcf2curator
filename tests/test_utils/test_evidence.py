"""Tests for supporting-read collection, SAME_START_READS and evidence BAMs."""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from vcf2curator.exceptions import InputValidationError
from vcf2curator.utils.evidence import collect_alignments, has_same_start_reads, write_evidence_bam


def _names(event: str) -> list[str]:
    return [f"{event}_r{i}" for i in range(1, 5)]


def test_collect_skips_supplementary_and_unwanted(sv_bam: Path) -> None:
    found = collect_alignments(sv_bam, set(_names("del")) | {"not_in_bam"})
    assert sorted(found) == _names("del")
    assert all(len(alns) == 2 for alns in found.values())
    assert all(not a.is_supplementary for alns in found.values() for a in alns)


def test_collect_bad_bam(fixtures_dir: Path) -> None:
    with pytest.raises(InputValidationError):
        collect_alignments(fixtures_dir / "garbage.txt", {"x"})


def test_same_start_reads(sv_bam: Path) -> None:
    found = collect_alignments(sv_bam, set(_names("del") + _names("dup")))
    assert not has_same_start_reads(_names("del"), found, fraction=0.6, min_distinct=3)
    assert has_same_start_reads(_names("dup"), found, fraction=0.6, min_distinct=3)
    # Too few distinct signatures fails even when none dominates.
    assert has_same_start_reads(_names("del"), found, fraction=0.6, min_distinct=5)
    # Reads missing from the BAM still count toward the cutoff denominator.
    assert has_same_start_reads(["del_r1", "ghost"], found, fraction=0.5, min_distinct=1)


def test_write_evidence_bam(sv_bam: Path, tmp_path: Path) -> None:
    found = collect_alignments(sv_bam, set(_names("del") + _names("inv")))
    out = tmp_path / "evidence" / "sample.evidence.bam"
    count = write_evidence_bam(sv_bam, found, _names("del"), out)
    assert count == 8
    assert Path(str(out) + ".bai").exists()
    assert not (out.parent / (out.name + ".unsorted.bam")).exists()
    with pysam.AlignmentFile(str(out), "rb") as bam:
        starts = [r.reference_start for r in bam.fetch("1")]
        names = {r.query_name for r in bam.fetch("1")}
    assert starts == sorted(starts)
    assert names == set(_names("del"))
