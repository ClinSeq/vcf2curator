"""Tests for the GRIDSS SV -> Curator table conversion."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from vcf2curator.annotations.cgc import load_cgc_by_symbol
from vcf2curator.annotations.curation import load_curation
from vcf2curator.annotations.intervals import load_bed
from vcf2curator.commands.svs import (
    SvAnnotations,
    SvFilterConfig,
    build_rows,
    drop_same_start_events,
    read_events,
)
from vcf2curator.curator.schema import SV_COLUMNS, VcfType
from vcf2curator.exceptions import InputValidationError
from vcf2curator.utils.evidence import collect_alignments


def _annotations(fixtures_dir: Path, *, panel: bool) -> SvAnnotations:
    return SvAnnotations(
        genes=load_bed(fixtures_dir / "sv_genes.bed", require_name=True),
        exons=load_bed(fixtures_dir / "sv_exons.bed"),
        cgc=load_cgc_by_symbol(fixtures_dir / "cgc.tsv"),
        targets=load_bed(fixtures_dir / "sv_targets.bed") if panel else None,
        curation=load_curation(fixtures_dir / "sv_curation.csv"),
    )


def _rows(fixtures_dir: Path, vcf: str, vcftype: VcfType, *, panel: bool) -> list[dict[str, Any]]:
    events = read_events(fixtures_dir / vcf, vcftype, need_read_names=False)
    rows = build_rows(
        events,
        _annotations(fixtures_dir, panel=panel),
        vcftype=vcftype,
        sdid="P-0001-T-01",
        project="PROJ",
        filters=SvFilterConfig(),
    )
    return [dict(r) for r in rows]


def test_read_events_skips_filtered_single_decoy_and_mates(fixtures_dir: Path) -> None:
    events = read_events(fixtures_dir / "somatic.gridss.vcf", VcfType.SOMATIC, False)
    # One record per event, both sides for the translocation.
    assert [(e.bnd.chrom, e.bnd.pos, e.svtype.value) for e in events] == [
        ("1", 1000, "DEL"),
        ("1", 3000, "DUP"),
        ("1", 5000, "INV"),
        ("1", 7000, "INS"),
        ("1", 8000, "TRA"),
        ("7", 500, "TRA"),
    ]


def test_somatic_panel_rows(fixtures_dir: Path) -> None:
    rows = _rows(fixtures_dir, "somatic.gridss.vcf", VcfType.SOMATIC, panel=True)
    assert all(tuple(r) == SV_COLUMNS for r in rows)
    # Legacy sort order: gene pair, then coordinates.
    assert [(r["SVTYPE"], r["CHROM_A"]) for r in rows] == [
        ("TRA", "1"),
        ("TRA", "7"),
        ("DEL", "1"),
        ("DUP", "1"),
        ("INV", "1"),
        ("INS", "1"),
    ]
    by_type = {r["SVTYPE"]: r for r in rows}

    deletion = by_type["DEL"]
    assert deletion == {
        "CHROM_A": "1",
        "START_A": 999,
        "END_A": 1000,
        "CHROM_B": "1",
        "START_B": 1999,
        "END_B": 2000,
        "IGV_COORD": "1:1000 1:1999",
        "SVTYPE": "DEL",
        "SV_LENGTH": 1001,
        "SUPPORT_READS": 30,
        "VAF": "0.25",
        "SDID": "P-0001-T-01",
        "SAMPLE": "somatic",
        "GENE_A": "BRCA1",
        "GENE_B": "GENEX",
        "GENE_A-GENE_B-sorted": "BRCA1,GENEX",
        "CGC_ANN": "TSG",
        "CURATOR": "YES",
        "CURATE": "BRCA1 check LOH",
    }
    assert by_type["INS"]["SV_LENGTH"] == 20
    # On target but no exon in the span -> NO.
    assert by_type["INV"]["CURATOR"] == "NO"
    assert by_type["DUP"]["CURATOR"] == "NO"

    tra = rows[0]
    assert (tra["START_A"], tra["END_A"], tra["START_B"], tra["END_B"]) == (8000, 8001, 499, 500)
    assert tra["SV_LENGTH"] == "NA"
    assert (tra["GENE_A"], tra["GENE_B"], tra["CGC_ANN"]) == ("NA", "BRAF", "oncogene")
    assert [r["CURATOR"] for r in rows[:2]] == ["YES", "YES"]


def test_somatic_wgs_uses_cgc(fixtures_dir: Path) -> None:
    rows = _rows(fixtures_dir, "somatic.gridss.vcf", VcfType.SOMATIC, panel=False)
    curated = {(r["SVTYPE"], r["CHROM_A"]) for r in rows if r["CURATOR"] == "YES"}
    assert curated == {("TRA", "1"), ("TRA", "7"), ("DEL", "1")}


def test_germline_rules_and_chr_prefix(fixtures_dir: Path) -> None:
    rows = _rows(fixtures_dir, "germline.gridss.vcf", VcfType.GERMLINE, panel=False)
    by_chrom = {r["CHROM_A"]: r for r in rows}
    assert sorted(by_chrom) == ["1", "7", "chr17"]
    assert by_chrom["1"]["CURATOR"] == "YES"
    assert by_chrom["1"]["VAF"] == "0.5"
    assert by_chrom["1"]["SAMPLE"] == "germline"
    # BRAF is in CGC with an exon hit, but only 10 supporting reads (< 40).
    assert by_chrom["7"]["CGC_ANN"] == "oncogene"
    assert by_chrom["7"]["CURATOR"] == "NO"
    # Germline curation entries do not match the somatic-only curation row.
    assert by_chrom["1"]["CURATE"] == ""


def test_drop_same_start_events(fixtures_dir: Path, sv_bam: Path) -> None:
    events = read_events(fixtures_dir / "somatic.gridss.vcf", VcfType.SOMATIC, True)
    alignments = collect_alignments(sv_bam, {n for e in events for n in e.read_names})
    kept = drop_same_start_events(events, alignments, SvFilterConfig())
    assert [e.svtype.value for e in kept] == ["DEL", "INV", "INS", "TRA", "TRA"]


def test_not_a_gridss_vcf(fixtures_dir: Path) -> None:
    with pytest.raises(InputValidationError, match="GRIDSS"):
        read_events(fixtures_dir / "somatic.vcf", VcfType.SOMATIC, False)


def test_missing_bpnames_header_when_bam_used(fixtures_dir: Path, tmp_path: Path) -> None:
    text = (fixtures_dir / "somatic.gridss.vcf").read_text()
    vcf = tmp_path / "no_bpnames.vcf"
    vcf.write_text("".join(line for line in text.splitlines(True) if "ID=BPNAMES" not in line))
    with pytest.raises(InputValidationError, match="BPNAMES"):
        read_events(vcf, VcfType.SOMATIC, True)


def test_unreadable_vcf(fixtures_dir: Path) -> None:
    with pytest.raises(InputValidationError):
        read_events(fixtures_dir / "garbage.txt", VcfType.SOMATIC, False)
