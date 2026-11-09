"""End-to-end CLI tests using Typer's CliRunner."""

from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

from vcf2curator.cli import app
from vcf2curator.curator.schema import (
    QC_PANEL_COLUMNS,
    QC_WGS_COLUMNS,
    SOMATIC_COLUMNS,
    SV_COLUMNS,
)

runner = CliRunner()


def test_version() -> None:
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert "vcf2curator" in result.stdout


def test_snvs_end_to_end(fixtures_dir: Path, tmp_path: Path) -> None:
    out = tmp_path / "out.curator.tsv"
    result = runner.invoke(
        app,
        [
            "snvs",
            "--input",
            str(fixtures_dir / "somatic.vcf"),
            "--output",
            str(out),
            "--vcftype",
            "somatic",
            "--oncokb",
            str(fixtures_dir / "oncokb.tsv"),
            "--cgc",
            str(fixtures_dir / "cgc.tsv"),
        ],
    )
    assert result.exit_code == 0, result.stdout
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "\t".join(SOMATIC_COLUMNS)
    assert len(lines) == 2  # header + one kept variant


def test_snvs_dry_run_writes_nothing(fixtures_dir: Path, tmp_path: Path) -> None:
    out = tmp_path / "out.curator.tsv"
    result = runner.invoke(
        app,
        [
            "snvs",
            "--input",
            str(fixtures_dir / "somatic.vcf"),
            "--output",
            str(out),
            "--vcftype",
            "somatic",
            "--oncokb",
            str(fixtures_dir / "oncokb.tsv"),
            "--dry-run",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert "dry-run" in result.stdout
    assert not out.exists()


def test_validate_ok(fixtures_dir: Path) -> None:
    result = runner.invoke(
        app,
        ["validate", "--input", str(fixtures_dir / "somatic.vcf")],
    )
    assert result.exit_code == 0


def test_validate_bad_input_exit_code(fixtures_dir: Path) -> None:
    result = runner.invoke(
        app,
        ["validate", "--input", str(fixtures_dir / "no_csq.vcf")],
    )
    assert result.exit_code == 3


def _svs_args(fixtures_dir: Path, out: Path, *extra: str) -> list[str]:
    return [
        "svs",
        "--input",
        str(fixtures_dir / "somatic.gridss.vcf"),
        "--output",
        str(out),
        "--vcftype",
        "somatic",
        "--sdid",
        "P-0001-T-01",
        "--project",
        "PROJ",
        "--genes",
        str(fixtures_dir / "sv_genes.bed"),
        "--exons",
        str(fixtures_dir / "sv_exons.bed"),
        "--cgc",
        str(fixtures_dir / "cgc.tsv"),
        *extra,
    ]


def test_svs_end_to_end_with_evidence_bam(fixtures_dir: Path, sv_bam: Path, tmp_path: Path) -> None:
    out = tmp_path / "out.sv.tsv"
    evidence = tmp_path / "out.evidence.bam"
    args = _svs_args(
        fixtures_dir,
        out,
        "--targets",
        str(fixtures_dir / "sv_targets.bed"),
        "--curation-ann",
        str(fixtures_dir / "sv_curation.csv"),
        "--bam",
        str(sv_bam),
        "--evidence-bam",
        str(evidence),
    )
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.stdout
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "\t".join(SV_COLUMNS)
    assert len(lines) == 1 + 5  # the DUP event fails SAME_START_READS
    assert evidence.exists()
    assert Path(str(evidence) + ".bai").exists()


def test_svs_dry_run_writes_nothing(fixtures_dir: Path, tmp_path: Path) -> None:
    out = tmp_path / "out.sv.tsv"
    result = runner.invoke(app, _svs_args(fixtures_dir, out, "--wgs", "--dry-run"))
    assert result.exit_code == 0, result.stdout
    assert "[dry-run] 6 rows" in result.stdout
    assert not out.exists()


@pytest.mark.parametrize(
    "extra",
    [
        [],  # neither --targets nor --wgs
        ["--wgs", "--targets", "TARGETS"],  # both
        ["--wgs", "--evidence-bam", "out.bam"],  # evidence BAM without --bam
    ],
)
def test_svs_bad_option_combinations(fixtures_dir: Path, tmp_path: Path, extra: list[str]) -> None:
    extra = [str(fixtures_dir / "sv_targets.bed") if e == "TARGETS" else e for e in extra]
    out = tmp_path / "out.sv.tsv"
    result = runner.invoke(app, _svs_args(fixtures_dir, out, *extra))
    assert result.exit_code == 3
    assert not out.exists()


def _qc_args(qc_dir: Path, role: str, metrics: str = "hsmetrics") -> list[str]:
    """Return the options for one sample, using the tumor/normal fixture files."""
    args = [
        f"--{role}-name",
        f"P-0001-{role[0].upper()}-01",
        f"--{role}-{metrics}",
        str(qc_dir / f"{role}.{metrics}.txt"),
        f"--{role}-insertsize",
        str(qc_dir / f"{role}.insertsize.txt"),
        f"--{role}-contest",
        str(qc_dir / f"{role}.contest.txt"),
    ]
    if metrics == "hsmetrics":
        args += [f"--{role}-markdup", str(qc_dir / f"{role}.markdup.txt")]
    return args


def _without(args: list[str], option: str) -> list[str]:
    i = args.index(option)
    return args[:i] + args[i + 2 :]


def test_qc_panel_tumor_normal(fixtures_dir: Path, tmp_path: Path) -> None:
    qc_dir = fixtures_dir / "qc"
    out = tmp_path / "case.qc.txt"
    args = [
        *_qc_args(qc_dir, "tumor"),
        "--tumor-msings",
        str(qc_dir / "tumor.MSI_Analysis.txt"),
        *_qc_args(qc_dir, "normal"),
    ]
    result = runner.invoke(app, ["qc", "-o", str(out), *args])
    assert result.exit_code == 0, result.stdout
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "\t".join(QC_PANEL_COLUMNS)
    assert [line.split("\t")[0] for line in lines[1:]] == ["P-0001-T-01", "P-0001-N-01"]
    assert lines[2].endswith("\tNA")  # normal has no mSINGS file


def test_qc_wgs_tumor_only(fixtures_dir: Path, tmp_path: Path) -> None:
    out = tmp_path / "case.qc.txt"
    args = _qc_args(fixtures_dir / "qc", "tumor", "wgsmetrics")
    result = runner.invoke(app, ["qc", "-o", str(out), *args])
    assert result.exit_code == 0, result.stdout
    assert out.read_text(encoding="utf-8").splitlines() == [
        "\t".join(QC_WGS_COLUMNS),
        "P-0001-T-01\t61.234\t1.21\t0.3\t168",
    ]


def test_qc_dry_run_writes_nothing(fixtures_dir: Path, tmp_path: Path) -> None:
    out = tmp_path / "case.qc.txt"
    args = ["qc", "-o", str(out), "--dry-run", *_qc_args(fixtures_dir / "qc", "tumor")]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.stdout
    assert "[dry-run] 1 rows" in result.stdout
    assert not out.exists()


@pytest.mark.parametrize(
    "case",
    [
        "no_tumor_name",
        "no_metrics",
        "both_metrics",
        "panel_without_markdup",
        "normal_incomplete",
        "mixed_modes",
        "markdup_with_wgs",
    ],
)
def test_qc_bad_option_combinations(fixtures_dir: Path, tmp_path: Path, case: str) -> None:
    qc_dir = fixtures_dir / "qc"
    tumor = _qc_args(qc_dir, "tumor")
    tumor_wgs = _qc_args(qc_dir, "tumor", "wgsmetrics")
    args = {
        "no_tumor_name": _without(tumor, "--tumor-name"),
        "no_metrics": _without(tumor, "--tumor-hsmetrics"),
        "both_metrics": [*tumor, "--tumor-wgsmetrics", str(qc_dir / "tumor.wgsmetrics.txt")],
        "panel_without_markdup": _without(tumor, "--tumor-markdup"),
        "normal_incomplete": [*tumor, "--normal-name", "P-0001-N-01"],
        "mixed_modes": [*tumor_wgs, *_qc_args(qc_dir, "normal")],
        "markdup_with_wgs": [*tumor_wgs, "--tumor-markdup", str(qc_dir / "tumor.markdup.txt")],
    }[case]
    out = tmp_path / "case.qc.txt"
    result = runner.invoke(app, ["qc", "-o", str(out), *args])
    assert result.exit_code == 3, result.stdout
    assert not out.exists()
