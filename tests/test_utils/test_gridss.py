"""Tests for GRIDSS breakend parsing and simple SV type classification."""

from __future__ import annotations

import pytest

from vcf2curator.utils.gridss import SvType, parse_breakend, simple_type, sv_length


@pytest.mark.parametrize(
    ("chrom", "pos", "ref", "alt", "expected"),
    [
        ("1", 1000, "A", "A[1:2000[", SvType.DEL),
        ("1", 2000, "C", "]1:1000]C", SvType.DEL),
        ("1", 3000, "A", "]1:4000]A", SvType.DUP),
        ("1", 4000, "C", "C[1:3000[", SvType.DUP),
        ("1", 5000, "A", "A]1:6000]", SvType.INV),
        ("1", 6000, "C", "[1:5000[C", SvType.INV),
        ("1", 7000, "A", "A" + "G" * 20 + "[1:7001[", SvType.INS),
        ("1", 7001, "C", "]1:7000]" + "G" * 20 + "C", SvType.INS),
        ("1", 8000, "A", "A[7:500[", SvType.TRA),
        ("chr1", 8000, "A", "A[chr7:500[", SvType.TRA),
    ],
)
def test_simple_type(chrom: str, pos: int, ref: str, alt: str, expected: SvType) -> None:
    bnd = parse_breakend(chrom, pos, ref, alt)
    assert bnd is not None
    assert simple_type(bnd) is expected


def test_parse_breakend_fields() -> None:
    bnd = parse_breakend("1", 1000, "A", "ATT[chr7:500[")
    assert bnd is not None
    assert (bnd.strand, bnd.partner_chrom, bnd.partner_pos, bnd.partner_strand) == (
        "+",
        "chr7",
        500,
        "-",
    )
    assert bnd.inserted == "TT"


@pytest.mark.parametrize("alt", ["A.", ".A", "<DEL>", "T"])
def test_single_or_non_breakend_alt_is_none(alt: str) -> None:
    assert parse_breakend("1", 100, "A", alt) is None


def test_deletion_with_small_insertion_stays_del() -> None:
    # 999 deleted bases, 5 inserted: far below the 0.7 insertion threshold.
    bnd = parse_breakend("1", 1000, "A", "AGGGGG[1:2000[")
    assert bnd is not None
    assert sv_length(bnd) == -999 + 5
    assert simple_type(bnd) is SvType.DEL


def test_sv_length_none_across_chromosomes() -> None:
    bnd = parse_breakend("1", 1000, "A", "A[7:500[")
    assert bnd is not None
    assert sv_length(bnd) is None
