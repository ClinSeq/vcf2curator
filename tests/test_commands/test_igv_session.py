"""Tests for IGV session XML generation."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from importlib.resources import files
from pathlib import Path

import pytest

from vcf2curator.commands.igv_session import (
    SNV_LAYOUT,
    SV_LAYOUT,
    PathStyle,
    SessionInputs,
    SessionLayout,
    bam_color,
    drop_empty_inputs,
    format_path,
    render_session,
)
from vcf2curator.exceptions import InputValidationError


def _touch(path: Path, content: str = "x") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path


@pytest.fixture
def inputs(tmp_path: Path) -> SessionInputs:
    data = tmp_path / "results"
    return SessionInputs(
        tumor_bam=_touch(data / "P-1-T-01.bam"),
        normal_bam=_touch(data / "P-1-N-01.bam"),
        tumor_evidence_bam=_touch(data / "P-1-T-01.evidence.bam"),
        normal_evidence_bam=_touch(data / "P-1-N-01.evidence.bam"),
        tumor_segments=_touch(data / "t_segments.bedGraph"),
        tumor_bins=_touch(data / "t_bins.bedGraph"),
        normal_segments=_touch(data / "n_segments.bedGraph"),
        normal_bins=_touch(data / "n_bins.bedGraph"),
        germline_afs=_touch(data / "afs.bedGraph"),
        somatic_vcf=_touch(data / "somatic.vep.vcf"),
        germline_vcf=_touch(data / "germline.vep.vcf"),
        targets=_touch(data / "targets.bed"),
    )


def _render(
    layout: SessionLayout,
    inputs: SessionInputs,
    outdir: Path,
    style: PathStyle = PathStyle.ABSOLUTE,
) -> ET.Element:
    template = files("vcf2curator.resources.igv").joinpath(layout.template).read_text()
    xml = render_session(
        template,
        layout,
        inputs,
        outdir=outdir,
        path_style=style,
        genome="hg38",
        reference_dir=Path("/ref/igv"),
    )
    return ET.fromstring(xml)


def _resources(root: ET.Element) -> list[str]:
    return [r.attrib["path"] for r in root.iter("Resource")]


def _data_track_names(root: ET.Element) -> list[str]:
    panel = next(p for p in root.iter("Panel") if p.get("name") == "DataPanel")
    return [t.attrib["name"] for t in panel.iter("Track")]


def _alignment_panels(root: ET.Element) -> list[str]:
    return [p.attrib["name"] for p in root.iter("Panel") if p.get("name", "").endswith(".bam")]


def test_snv_session_tracks_in_order(inputs: SessionInputs, tmp_path: Path) -> None:
    root = _render(SNV_LAYOUT, inputs, tmp_path / "out")
    assert _alignment_panels(root) == ["P-1-T-01.bam", "P-1-N-01.bam"]
    assert _data_track_names(root) == [
        "gDNA_segments",
        "gDNA_bins",
        "tDNA_segments",
        "tDNA_bins",
        "SNP_tDNA",
        "GERMLINE VEP VCF",
        "SOMATIC VEP VCF",
    ]


def test_sv_session_has_evidence_bams_and_no_vcfs(inputs: SessionInputs, tmp_path: Path) -> None:
    root = _render(SV_LAYOUT, inputs, tmp_path / "out")
    assert _alignment_panels(root) == [
        "P-1-T-01.bam",
        "P-1-N-01.bam",
        "P-1-N-01.evidence.bam",
        "P-1-T-01.evidence.bam",
    ]
    assert "SOMATIC VEP VCF" not in _data_track_names(root)
    assert not any(p.endswith(".vcf") for p in _resources(root))


def test_absolute_paths_and_options(inputs: SessionInputs, tmp_path: Path) -> None:
    root = _render(SNV_LAYOUT, inputs, tmp_path / "out")
    resources = _resources(root)
    assert resources[0] == str(inputs.targets)  # capture BED first, as legacy
    assert str(inputs.tumor_bam) in resources
    assert "/ref/igv/consensusBlacklist.bed" in resources
    assert root.attrib["genome"] == "hg38"
    track_ids = [t.attrib["id"] for t in root.iter("Track")]
    assert "hg38_genes" in track_ids
    assert f"{inputs.tumor_bam}_coverage" in track_ids


def test_relative_paths(inputs: SessionInputs, tmp_path: Path) -> None:
    root = _render(SNV_LAYOUT, inputs, tmp_path / "out", PathStyle.RELATIVE)
    assert "../results/P-1-T-01.bam" in _resources(root)
    # Reference tracks stay absolute.
    assert "/ref/igv/simple_repeats.bed" in _resources(root)


def test_missing_roles_are_left_out(tmp_path: Path) -> None:
    only_bam = SessionInputs(tumor_bam=_touch(tmp_path / "t.bam"))
    root = _render(SNV_LAYOUT, only_bam, tmp_path)
    assert _alignment_panels(root) == ["t.bam"]
    assert _data_track_names(root) == []


def test_paths_are_xml_escaped(tmp_path: Path) -> None:
    bam = _touch(tmp_path / 'a&b"c.bam')
    root = _render(SV_LAYOUT, SessionInputs(tumor_bam=bam), tmp_path)
    assert str(bam) in _resources(root)


def test_drop_empty_inputs_warns(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    empty = _touch(tmp_path / "empty.bedGraph", "")
    bam = _touch(tmp_path / "t.bam")
    kept = drop_empty_inputs(SessionInputs(tumor_bam=bam, tumor_bins=empty))
    assert kept == SessionInputs(tumor_bam=bam)
    assert "--tumor-bins" in caplog.text


def test_invalid_template_placeholder(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="placeholder"):
        render_session(
            "<Session>{unknown}</Session>",
            SNV_LAYOUT,
            SessionInputs(),
            outdir=tmp_path,
            path_style=PathStyle.ABSOLUTE,
            genome="hg19",
            reference_dir=tmp_path,
        )


@pytest.mark.parametrize(
    ("name", "color"),
    [
        ("P-1-T-01_DEL.bam", "53,116,199"),
        ("P-1-CFDNA-01_TRA.bam", "239,133,62"),
        ("P-1-N-01_INV.bam", "69,136,51"),
        ("P-1-T-01.bam", "175,175,175"),
    ],
)
def test_bam_color(name: str, color: str) -> None:
    assert bam_color(name) == color


def test_format_path_relative_to_sibling(tmp_path: Path) -> None:
    file = _touch(tmp_path / "a" / "x.bam")
    assert format_path(file, tmp_path / "b", PathStyle.RELATIVE) == "../a/x.bam"
