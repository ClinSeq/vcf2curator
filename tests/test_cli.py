"""End-to-end CLI tests using Typer's CliRunner."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from vcf2curator.cli import app
from vcf2curator.curator.schema import SOMATIC_COLUMNS

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
