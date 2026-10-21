"""Tests for the SNV/indel VCF parser."""

from __future__ import annotations

from pathlib import Path

import pytest

from vcf2curator.annotations.aminoacid import tri_to_single_aa
from vcf2curator.commands.snvs import SnvParser
from vcf2curator.curator.schema import GERMLINE_COLUMNS, SOMATIC_COLUMNS, VcfType
from vcf2curator.exceptions import InputValidationError


def _somatic_parser(fixtures_dir: Path) -> SnvParser:
    return SnvParser(
        vcftype=VcfType.SOMATIC,
        oncokb_path=fixtures_dir / "oncokb.tsv",
        cgc_path=fixtures_dir / "cgc.tsv",
    )


def _germline_parser(fixtures_dir: Path) -> SnvParser:
    return SnvParser(
        vcftype=VcfType.GERMLINE,
        oncokb_path=fixtures_dir / "oncokb.tsv",
        cgc_path=fixtures_dir / "cgc.tsv",
    )


def test_tri_to_single_aa() -> None:
    assert tri_to_single_aa("Val600Glu") == "V600E"
    assert tri_to_single_aa("Gln41Ter") == "Q41Ter"


def test_somatic_happy_path(fixtures_dir: Path) -> None:
    rows = list(_somatic_parser(fixtures_dir).parse(fixtures_dir / "somatic.vcf"))
    # Only the MODERATE-impact BRAF variant survives filtering; the LOW one drops.
    assert len(rows) == 1
    row = rows[0]
    assert set(row) == set(SOMATIC_COLUMNS)
    assert row["CHROM"] == "7"
    assert row["START"] == 140453135
    assert row["END"] == 140453136
    assert row["REF"] == "A"
    assert row["ALT"] == "T"
    assert row["GENE"] == "BRAF"
    assert row["IMPACT"] == "MODERATE"
    assert row["AMINO_ACID"] == "V600E"
    assert row["T_DP"] == 120
    assert row["T_ALT"] == 40
    assert row["T_VAF"] == pytest.approx(0.33)
    assert row["N_DP"] == 60
    assert row["N_ALT"] == 2
    assert row["CGC_ANN"] == "oncogene"
    assert row["NUM_TOOLS"] == 2
    assert row["RSID"] == "rs113488022"
    # OncoKB is emitted as the raw [oncogenicity, effect, pmids] entry.
    assert row["OncoKB"] == ["Oncogenic", "Gain-of-function", "12345"]


def test_germline_happy_path(fixtures_dir: Path) -> None:
    rows = list(_germline_parser(fixtures_dir).parse(fixtures_dir / "germline.vcf"))
    # HIGH-impact stop_gained survives; non-pathogenic missense is filtered out.
    assert len(rows) == 1
    row = rows[0]
    assert set(row) == set(GERMLINE_COLUMNS)
    assert row["CHROM"] == "17"
    assert row["GENE"] == "BRCA1"
    assert row["IMPACT"] == "HIGH"
    assert row["N_DP"] == 30
    assert row["N_ALT"] == 15
    assert row["N_VAF"] == pytest.approx(0.5)
    assert row["T_VAF"] == pytest.approx(0.0)
    assert row["RSID"] == "rs80357065"
    assert row["CGC_ANN"] == "TSG"


def test_parse_without_csq_raises(fixtures_dir: Path) -> None:
    parser = _somatic_parser(fixtures_dir)
    with pytest.raises(InputValidationError, match="CSQ"):
        list(parser.parse(fixtures_dir / "no_csq.vcf"))


def test_open_garbage_raises(fixtures_dir: Path) -> None:
    parser = _somatic_parser(fixtures_dir)
    with pytest.raises(InputValidationError):
        list(parser.parse(fixtures_dir / "garbage.txt"))


def test_validate_ok(fixtures_dir: Path) -> None:
    parser = _somatic_parser(fixtures_dir)
    parser.validate(fixtures_dir / "somatic.vcf")  # should not raise


def test_validate_rejects_missing_csq(fixtures_dir: Path) -> None:
    parser = _somatic_parser(fixtures_dir)
    with pytest.raises(InputValidationError):
        parser.validate(fixtures_dir / "no_csq.vcf")
