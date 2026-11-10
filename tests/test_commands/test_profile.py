"""Tests for the profile command."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest
import yaml
from typer.testing import CliRunner

from vcf2curator.cli import app
from vcf2curator.commands.profile import (
    PathStyle,
    ProfileInputs,
    build_profile,
    check_qc_samples,
    check_table_header,
    format_path,
)
from vcf2curator.curator.schema import (
    GERMLINE_COLUMNS,
    QC_PANEL_COLUMNS,
    QC_WGS_COLUMNS,
    SOMATIC_COLUMNS,
    SV_COLUMNS,
)
from vcf2curator.exceptions import InputValidationError

runner = CliRunner()


def _table(path: Path, columns: tuple[str, ...], *rows: dict[str, str]) -> Path:
    lines = ["\t".join(columns)]
    lines += ["\t".join(row.get(col, "") for col in columns) for row in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def _qc(path: Path, columns: tuple[str, ...] = QC_PANEL_COLUMNS) -> Path:
    return _table(path, columns, {"SAMP": "P-0001-T-01"}, {"SAMP": "P-0001-N-01"})


@dataclass
class Case:
    """A complete set of synthetic case files."""

    root: Path
    snvs_somatic: Path
    snvs_germline: Path
    svs_somatic: Path
    svs_germline: Path
    evidence_bam: Path
    qc: Path
    multiqc: Path
    session_snvs: Path
    session_svs: Path
    snapshots_snvs: Path
    snapshots_svs: Path
    frankenplot_html: Path
    frankenplot_images: tuple[Path, Path]


def _touch(path: Path, text: str = "x") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


@pytest.fixture
def case(tmp_path: Path) -> Case:
    root = tmp_path / "results"
    root.mkdir()
    snap_snvs = root / "snapshots" / "snvs"
    snap_svs = root / "snapshots" / "svs"
    _touch(snap_snvs / "chr7_1-2.png")
    snap_svs.mkdir(parents=True)
    bam = _touch(root / "evidence.bam")
    _touch(root / "evidence.bam.bai")
    return Case(
        root=root,
        snvs_somatic=_table(root / "somatic.tsv", SOMATIC_COLUMNS),
        snvs_germline=_table(root / "germline.tsv", GERMLINE_COLUMNS),
        svs_somatic=_table(root / "svs.somatic.tsv", SV_COLUMNS),
        svs_germline=_table(root / "svs.germline.tsv", SV_COLUMNS),
        evidence_bam=bam,
        qc=_qc(root / "qc.tsv"),
        multiqc=_touch(root / "multiqc_report.html"),
        session_snvs=_touch(root / "igv_session_snps.xml"),
        session_svs=_touch(root / "igv_session_sv.xml"),
        snapshots_snvs=snap_snvs,
        snapshots_svs=snap_svs,
        frankenplot_html=_touch(root / "frankenplot.html"),
        frankenplot_images=(_touch(root / "fp1.png"), _touch(root / "fp2.png")),
    )


def _required(case: Case, output: Path) -> list[str]:
    return [
        "profile",
        "--case-id",
        "P-0001",
        "--tumor-name",
        "P-0001-T-01",
        "--snvs-somatic",
        str(case.snvs_somatic),
        "--qc",
        str(case.qc),
        "-o",
        str(output),
    ]


def _full(case: Case, output: Path) -> list[str]:
    return [
        *_required(case, output),
        "--project",
        "PROBIO",
        "--cancer-type",
        "prostate",
        "--normal-name",
        "P-0001-N-01",
        "--snvs-germline",
        str(case.snvs_germline),
        "--svs-somatic",
        str(case.svs_somatic),
        "--svs-germline",
        str(case.svs_germline),
        "--sv-evidence-bam",
        str(case.evidence_bam),
        "--multiqc",
        str(case.multiqc),
        "--igv-session-snvs",
        str(case.session_snvs),
        "--igv-session-svs",
        str(case.session_svs),
        "--igv-snapshots-snvs",
        str(case.snapshots_snvs),
        "--igv-snapshots-svs",
        str(case.snapshots_svs),
        "--frankenplot-html",
        str(case.frankenplot_html),
        "--frankenplot-image",
        str(case.frankenplot_images[0]),
        "--frankenplot-image",
        str(case.frankenplot_images[1]),
    ]


def _load(path: Path) -> Any:
    return yaml.safe_load(path.read_text())


# -- unit -------------------------------------------------------------------


def test_check_table_header_accepts_qc_variants(tmp_path: Path) -> None:
    allowed = (QC_PANEL_COLUMNS, QC_WGS_COLUMNS)
    assert check_table_header(_qc(tmp_path / "p.tsv"), allowed) == QC_PANEL_COLUMNS
    assert check_table_header(_qc(tmp_path / "w.tsv", QC_WGS_COLUMNS), allowed) == QC_WGS_COLUMNS


def test_check_table_header_rejects_wrong_schema(tmp_path: Path) -> None:
    germline = _table(tmp_path / "g.tsv", GERMLINE_COLUMNS)
    with pytest.raises(InputValidationError, match="header does not match"):
        check_table_header(germline, (SOMATIC_COLUMNS,))


def test_check_table_header_rejects_empty(tmp_path: Path) -> None:
    empty = tmp_path / "e.tsv"
    empty.write_text("")
    with pytest.raises(InputValidationError, match="empty table"):
        check_table_header(empty, (SV_COLUMNS,))


def test_check_qc_samples(tmp_path: Path) -> None:
    qc = _qc(tmp_path / "qc.tsv")
    check_qc_samples(qc, ["P-0001-T-01", "P-0001-N-01"])
    with pytest.raises(InputValidationError, match="P-0002-T-01"):
        check_qc_samples(qc, ["P-0002-T-01"])


def test_format_path_relative_keeps_symlinks(tmp_path: Path) -> None:
    target = _touch(tmp_path / "work" / "real.tsv")
    link = tmp_path / "results" / "link.tsv"
    link.parent.mkdir()
    link.symlink_to(target)
    assert format_path(link, tmp_path / "results", PathStyle.RELATIVE) == "link.tsv"
    assert format_path(link, tmp_path, PathStyle.ABSOLUTE) == str(link)


def test_build_profile_omits_missing(case: Case) -> None:
    inputs = ProfileInputs(
        case_id="P-0001", tumor_name="T", snvs_somatic=case.snvs_somatic, qc=case.qc
    )
    assert build_profile(inputs, base=case.root, style=PathStyle.RELATIVE) == {
        "schema_version": 1,
        "case_id": "P-0001",
        "samples": {"tumor": "T"},
        "variants": {"snvs": {"somatic": "somatic.tsv"}},
        "qc": {"table": "qc.tsv"},
    }


# -- CLI --------------------------------------------------------------------


def test_cli_minimal(case: Case, tmp_path: Path) -> None:
    output = tmp_path / "out" / "profile.yaml"
    result = runner.invoke(app, _required(case, output))
    assert result.exit_code == 0, result.output
    assert _load(output) == {
        "schema_version": 1,
        "case_id": "P-0001",
        "samples": {"tumor": "P-0001-T-01"},
        "variants": {"snvs": {"somatic": str(case.snvs_somatic)}},
        "qc": {"table": str(case.qc)},
    }


def test_cli_full(case: Case, tmp_path: Path) -> None:
    output = tmp_path / "profile.yaml"
    result = runner.invoke(app, _full(case, output))
    assert result.exit_code == 0, result.output
    data = _load(output)
    assert list(data) == [
        "schema_version",
        "case_id",
        "project",
        "cancer_type",
        "samples",
        "variants",
        "qc",
        "igv",
        "cnv",
    ]
    assert data["samples"] == {"tumor": "P-0001-T-01", "normal": "P-0001-N-01"}
    assert data["variants"]["svs"]["evidence_bam"] == str(case.evidence_bam)
    assert data["igv"]["snapshots"]["svs"] == str(case.snapshots_svs)
    assert data["cnv"]["frankenplot_images"] == [str(p) for p in case.frankenplot_images]


def test_cli_relative_paths(case: Case) -> None:
    output = case.root / "curator" / "profile.yaml"
    result = runner.invoke(app, [*_full(case, output), "--path-style", "relative"])
    assert result.exit_code == 0, result.output
    data = _load(output)
    assert data["variants"]["snvs"]["somatic"] == "../somatic.tsv"
    assert data["igv"]["snapshots"]["snvs"] == "../snapshots/snvs"


def test_cli_warns_on_empty_snapshot_dir(
    case: Case, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    with caplog.at_level(logging.WARNING, logger="vcf2curator"):
        result = runner.invoke(app, _full(case, tmp_path / "profile.yaml"))
    assert result.exit_code == 0, result.output
    assert f"No PNG snapshots in {case.snapshots_svs}" in caplog.text


def test_cli_germline_as_somatic_fails(case: Case, tmp_path: Path) -> None:
    args = _required(case, tmp_path / "profile.yaml")
    args[args.index("--snvs-somatic") + 1] = str(case.snvs_germline)
    result = runner.invoke(app, args)
    assert result.exit_code == 3


def test_cli_bad_sv_header_fails(case: Case, tmp_path: Path) -> None:
    bad = _table(tmp_path / "bad_svs.tsv", SV_COLUMNS[:-1])
    args = [*_required(case, tmp_path / "profile.yaml"), "--svs-somatic", str(bad)]
    assert runner.invoke(app, args).exit_code == 3


def test_cli_qc_sample_mismatch_fails(case: Case, tmp_path: Path) -> None:
    args = [*_required(case, tmp_path / "profile.yaml"), "--normal-name", "OTHER-N"]
    assert runner.invoke(app, args).exit_code == 3


def test_cli_empty_file_fails(case: Case, tmp_path: Path) -> None:
    empty = _touch(tmp_path / "multiqc.html", "")
    args = [*_required(case, tmp_path / "profile.yaml"), "--multiqc", str(empty)]
    assert runner.invoke(app, args).exit_code == 3


def test_cli_dry_run_writes_nothing(case: Case, tmp_path: Path) -> None:
    output = tmp_path / "out" / "profile.yaml"
    result = runner.invoke(app, [*_required(case, output), "--dry-run"])
    assert result.exit_code == 0, result.output
    assert not output.parent.exists()


def test_cli_missing_required_option(case: Case, tmp_path: Path) -> None:
    args = _required(case, tmp_path / "profile.yaml")
    del args[args.index("--qc") : args.index("--qc") + 2]
    assert runner.invoke(app, args).exit_code == 2
