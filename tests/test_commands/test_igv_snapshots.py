"""Tests for the igv-snapshots command."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

from vcf2curator.cli import app
from vcf2curator.commands.igv_snapshots import (
    SNV_BATCH_NAME,
    SNV_REQUIRED_COLUMNS,
    SV_BATCH_NAME,
    batch_script,
    read_table,
    snv_commands,
    snv_goto_pos,
    sv_commands,
    ucsc_chrom,
)
from vcf2curator.exceptions import InputValidationError

runner = CliRunner()


@pytest.fixture
def igv_dir(fixtures_dir: Path) -> Path:
    return fixtures_dir / "igv"


class _FakeRun:
    """Records subprocess.run calls and returns a fixed exit code."""

    def __init__(self, returncode: int = 0) -> None:
        self.returncode = returncode
        self.calls: list[list[str]] = []

    def __call__(self, cmd: list[str], **_: Any) -> subprocess.CompletedProcess[bytes]:
        self.calls.append(cmd)
        return subprocess.CompletedProcess(cmd, self.returncode)


@pytest.fixture
def fake_run(monkeypatch: pytest.MonkeyPatch) -> _FakeRun:
    fake = _FakeRun()
    monkeypatch.setattr(subprocess, "run", fake)
    return fake


# -- pure functions ---------------------------------------------------------


@pytest.mark.parametrize(("chrom", "expected"), [("7", "chr7"), ("chr7", "chr7"), ("X", "chrX")])
def test_ucsc_chrom(chrom: str, expected: str) -> None:
    assert ucsc_chrom(chrom) == expected


@pytest.mark.parametrize(
    ("ref", "alt", "expected"),
    [
        ("A", "T", 100),  # SNV
        ("A", "AT", 100),  # insertion
        ("CTTG", "C", 102),  # deletion: END + round(4 / 2)
        ("CTT", "C", 102),  # deletion: END + round(3 / 2)
        ("G", "GA,GAA", 100),  # multi-allelic insertion
        ("AT", "A,ATT", 101),  # first ALT shorter -> deletion
    ],
)
def test_snv_goto_pos(ref: str, alt: str, expected: int) -> None:
    assert snv_goto_pos(ref, alt, 100) == expected


def test_snv_commands(igv_dir: Path) -> None:
    rows = read_table(igv_dir / "snvs.tsv", SNV_REQUIRED_COLUMNS)
    assert snv_commands(rows) == [
        "goto chr7:140453136",
        "sort base",
        "snapshot chr7_140453135-140453136.png",
        "goto chr17:7578411",
        "sort base",
        "snapshot chr17_7578405-7578409.png",
        "goto chr12:25398281",
        "sort base",
        "snapshot chr12_25398280-25398281.png",
    ]


def test_snv_commands_skips_duplicate_names(caplog: pytest.LogCaptureFixture) -> None:
    row = {"CHROM": "1", "START": "9", "END": "10", "REF": "A", "ALT": "T"}
    assert len(snv_commands([row, dict(row)])) == 3
    assert "Duplicate snapshot chr1_9-10.png" in caplog.text


def test_snv_commands_rejects_non_integer() -> None:
    row = {"CHROM": "1", "START": "x", "END": "10", "REF": "A", "ALT": "T"}
    with pytest.raises(InputValidationError, match="START must be an integer"):
        snv_commands([row])


def test_sv_commands_skips_curator_no() -> None:
    rows = [
        {
            "CHROM_A": "21",
            "START_A": "999",
            "END_A": "1000",
            "CHROM_B": "chr21",
            "START_B": "4999",
            "END_B": "5000",
            "CURATOR": "YES",
        },
        {
            "CHROM_A": "1",
            "START_A": "1",
            "END_A": "2",
            "CHROM_B": "2",
            "START_B": "3",
            "END_B": "4",
            "CURATOR": "NO",
        },
    ]
    assert sv_commands(rows) == [
        "squish tb_simpleRepeat",
        "goto chr21:850-1150",
        "sort base",
        "snapshot chr21_999-1000.png",
        "goto chr21:4850-5150",
        "sort base",
        "snapshot chr21_4999-5000.png",
    ]


def test_batch_script(tmp_path: Path) -> None:
    session = tmp_path / "s.xml"
    script = batch_script(["goto chr1:1"], genome="hg38", outdir=tmp_path, session=session)
    assert script.splitlines() == [
        "new",
        "genome hg38",
        f"snapshotDirectory {tmp_path.resolve()}",
        f"load {session.resolve()}",
        "goto chr1:1",
        "exit",
    ]


def test_read_table_missing_column(tmp_path: Path) -> None:
    table = tmp_path / "t.tsv"
    table.write_text("CHROM\tSTART\n1\t2\n")
    with pytest.raises(InputValidationError, match="missing columns"):
        read_table(table, SNV_REQUIRED_COLUMNS)


def test_read_table_empty_file(tmp_path: Path) -> None:
    table = tmp_path / "t.tsv"
    table.write_text("")
    with pytest.raises(InputValidationError, match="no header"):
        read_table(table, SNV_REQUIRED_COLUMNS)


# -- CLI --------------------------------------------------------------------


def _args(sub: str, table: Path, session: Path, outdir: Path, *extra: str) -> list[str]:
    return [
        "igv-snapshots",
        sub,
        "-i",
        str(table),
        "--session",
        str(session),
        "-o",
        str(outdir),
        "--port",
        "60151",
        *extra,
    ]


def test_cli_snvs(igv_dir: Path, tmp_path: Path, fake_run: _FakeRun) -> None:
    outdir = tmp_path / "out"
    result = runner.invoke(
        app, _args("snvs", igv_dir / "snvs.tsv", igv_dir / "session.xml", outdir)
    )
    assert result.exit_code == 0, result.output

    batch = outdir.resolve() / SNV_BATCH_NAME
    script = batch.read_text().splitlines()
    assert script[:4] == [
        "new",
        "genome hg19",
        f"snapshotDirectory {outdir.resolve()}",
        f"load {(igv_dir / 'session.xml').resolve()}",
    ]
    assert script[-1] == "exit"
    assert sum(line.startswith("snapshot ") for line in script) == 3

    (cmd,) = fake_run.calls
    assert cmd[:3] == ["igv.sh", "-b", str(batch)]
    assert cmd[3] == "--preferences"
    assert Path(cmd[4]).name == "igvsnapshot.prefs"
    assert cmd[5:] == ["-p", "60151", "--igvDirectory", str(outdir.resolve())]


def test_cli_svs(igv_dir: Path, tmp_path: Path, fake_run: _FakeRun) -> None:
    outdir = tmp_path / "out"
    prefs = tmp_path / "my.prefs"
    prefs.write_text("SAM.SHOW_SOFT_CLIPPED=true\n")
    result = runner.invoke(
        app,
        _args(
            "svs",
            igv_dir / "svs.tsv",
            igv_dir / "session.xml",
            outdir,
            "--preferences",
            str(prefs),
            "--igv",
            "/opt/igv/igv.sh",
        ),
    )
    assert result.exit_code == 0, result.output

    script = (outdir / SV_BATCH_NAME).read_text()
    assert "goto chr21:39956720-39957020" in script
    assert "snapshot chr21_42880007-42880008.png" in script
    assert "chr1_1000-1001" not in script  # CURATOR=NO

    (cmd,) = fake_run.calls
    assert cmd[0] == "/opt/igv/igv.sh"
    assert cmd[4] == str(prefs.resolve())


def test_cli_dry_run_writes_nothing(igv_dir: Path, tmp_path: Path, fake_run: _FakeRun) -> None:
    outdir = tmp_path / "out"
    result = runner.invoke(
        app,
        _args("snvs", igv_dir / "snvs.tsv", igv_dir / "session.xml", outdir, "--dry-run"),
    )
    assert result.exit_code == 0, result.output
    assert not outdir.exists()
    assert fake_run.calls == []


def test_cli_empty_table_skips_igv(igv_dir: Path, tmp_path: Path, fake_run: _FakeRun) -> None:
    table = tmp_path / "empty.tsv"
    table.write_text("\t".join(SNV_REQUIRED_COLUMNS) + "\n")
    outdir = tmp_path / "out"
    result = runner.invoke(app, _args("snvs", table, igv_dir / "session.xml", outdir))
    assert result.exit_code == 0, result.output
    assert (outdir / SNV_BATCH_NAME).exists()
    assert fake_run.calls == []


def test_cli_missing_column_exit_code(igv_dir: Path, tmp_path: Path, fake_run: _FakeRun) -> None:
    result = runner.invoke(
        app, _args("svs", igv_dir / "snvs.tsv", igv_dir / "session.xml", tmp_path / "out")
    )
    assert result.exit_code == 3
    assert fake_run.calls == []


def test_cli_igv_failure_exit_code(
    igv_dir: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(subprocess, "run", _FakeRun(returncode=1))
    result = runner.invoke(
        app, _args("snvs", igv_dir / "snvs.tsv", igv_dir / "session.xml", tmp_path / "out")
    )
    assert result.exit_code == 6


def test_cli_igv_not_found_exit_code(igv_dir: Path, tmp_path: Path) -> None:
    result = runner.invoke(
        app,
        _args(
            "snvs",
            igv_dir / "snvs.tsv",
            igv_dir / "session.xml",
            tmp_path / "out",
            "--igv",
            str(tmp_path / "no-such-igv.sh"),
        ),
    )
    assert result.exit_code == 6
