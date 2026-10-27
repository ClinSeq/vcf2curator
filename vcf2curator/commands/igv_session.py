"""``igv-session`` command: per-sample input files -> IGV session XML files.

Replaces the session part of the legacy ``generate_symlinks.py``. Writes
``igv_session_snps.xml`` and ``igv_session_sv.xml`` into ``--outdir``.

Behavioural notes vs. the legacy script (intentional, flagged for review):

* Inputs are given per role on the command line; no symlinks, no directory
  scanning and no filename-regex classification.
* Tracks are named by role (``tDNA_segments``, ``SOMATIC VEP VCF``, ...), not
  derived from file names; BAM tracks keep the file name.
* Only tumor/normal BAMs and GRIDSS evidence BAMs are supported. nodups,
  markdups, clipoverlap selection, ``.mut`` and ``.gtf`` tracks and the WGS
  mode are dropped. Omit ``--somatic-vcf``/``--germline-vcf`` to leave VCFs out.
* Genome and reference directory are options instead of fixed ``hg19`` and
  ``/nfs/PROBIO/for_igv``.
* Missing input files fail; empty input files are skipped with a warning.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field, fields
from enum import StrEnum
from importlib.resources import files
from pathlib import Path
from typing import Annotated, Any
from xml.sax.saxutils import escape

import typer

from vcf2curator.commands._shared import QuietOpt, VerboseOpt, fail
from vcf2curator.exceptions import InputValidationError, Vcf2CuratorError
from vcf2curator.logging_config import configure_logging, get_logger

logger = get_logger()

SNV_SESSION_NAME = "igv_session_snps.xml"
SV_SESSION_NAME = "igv_session_sv.xml"

_TEMPLATES = files("vcf2curator.resources.igv")
_SNV_TEMPLATE = "igv_session_snvs.xml"
_SV_TEMPLATE = "igv_session_svs.xml"

# Alignment track colour by SV type in the BAM file name; grey otherwise.
_BAM_SV_TYPE = re.compile(r".*-(N|CFDNA|T)-.*(DEL|DUP|INV|TRA)\.bam$")
_SV_TYPE_COLORS = {
    "DEL": "53,116,199",
    "TRA": "239,133,62",
    "DUP": "204,55,48",
    "INV": "69,136,51",
}
_DEFAULT_BAM_COLOR = "175,175,175"

_RESOURCE = '<Resource path="{path}"/>'

_PANEL = """<Panel height="{height}" width="{width}" name="{name}">
        {tracks}
    </Panel>"""

_ALIGNMENT_TRACKS = """<Track autoScale="true" clazz="org.broad.igv.sam.CoverageTrack" color="{color}" colorScale="ContinuousColorScale;0.0;927.0;255,255,255;{color}" fontSize="10" id="{path}_coverage" name="{name} Coverage" snpThreshold="0.0001" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="1503.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track clazz="org.broad.igv.sam.SpliceJunctionTrack" fontSize="10" height="60" id="{path}_junctions" name="{name} Junctions" visible="false"/>
        <Track clazz="org.broad.igv.sam.AlignmentTrack" displayMode="EXPANDED" experimentType="OTHER" fontSize="10" id="{path}" name="{name}" visible="true">
            <RenderOptions colorOption="READ_STRAND"/>
        </Track>"""  # noqa: E501 - IGV track XML kept on one line per element

_SEGMENTS_TRACK = """<Track autoScale="false" clazz="org.broad.igv.track.DataSourceTrack" displayMode="EXPANDED" fontSize="10" height="80" id="{path}" name="{name}" altColor="0,0,255" color="255,0,0" renderer="BAR_CHART" visible="true" windowFunction="mean">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="2.0" minimum="-1.0" type="LINEAR"/>
        </Track>"""  # noqa: E501

_BINS_TRACK = """<Track autoScale="false" clazz="org.broad.igv.track.DataSourceTrack" displayMode="EXPANDED" fontSize="10" height="80" id="{path}" name="{name}" altColor="0,0,255" color="255,0,0" renderer="SCATTER_PLOT" visible="true" windowFunction="mean">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="2.0" minimum="-1.0" type="LINEAR"/>
        </Track>"""  # noqa: E501

_AFS_TRACK = """<Track clazz="org.broad.igv.track.DataSourceTrack" displayMode="EXPANDED" fontSize="10" height="80" id="{path}" name="{name}" renderer="SCATTER_PLOT" color="0,0,0" visible="true" windowFunction="mean">
            <DataRange baseline="0.5" drawBaseline="true" flipAxis="false" maximum="1" minimum="0" type="LINEAR"/>
        </Track>"""  # noqa: E501

_VCF_TRACK = '<Track clazz="org.broad.igv.variant.VariantTrack" color="0,0,178" displayMode="EXPANDED" fontSize="10" id="{path}" name="{name}" siteColorMode="ALLELE_FREQUENCY" squishedHeight="1" visible="true"/>'  # noqa: E501

_CAPTURE_BED_TRACK = '<Track clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" colorScale="ContinuousColorScale;0.0;103.0;255,255,255;0,0,178" fontSize="10" id="{path}" name="{name}" visible="true"/>'  # noqa: E501


class PathStyle(StrEnum):
    """How input file paths are written into the session."""

    ABSOLUTE = "absolute"
    RELATIVE = "relative"


class TrackKind(StrEnum):
    """Track layout used for an input file."""

    ALIGNMENT = "alignment"
    SEGMENTS = "segments"
    BINS = "bins"
    AFS = "afs"
    VCF = "vcf"


_DATA_TRACK_TEMPLATES = {
    TrackKind.SEGMENTS: _SEGMENTS_TRACK,
    TrackKind.BINS: _BINS_TRACK,
    TrackKind.AFS: _AFS_TRACK,
    TrackKind.VCF: _VCF_TRACK,
}


@dataclass(frozen=True)
class SessionInputs:
    """Per-sample input files; ``None`` means not given (or skipped as empty).

    Attributes:
        tumor_bam: Tumor alignment BAM.
        normal_bam: Normal alignment BAM.
        tumor_evidence_bam: Tumor GRIDSS evidence BAM (SV session only).
        normal_evidence_bam: Normal GRIDSS evidence BAM (SV session only).
        tumor_segments: Tumor CNV segments bedGraph.
        tumor_bins: Tumor CNV bins (profile) bedGraph.
        normal_segments: Normal CNV segments bedGraph.
        normal_bins: Normal CNV bins (profile) bedGraph.
        germline_afs: Germline SNP allele frequencies in the tumor (bedGraph).
        somatic_vcf: Somatic VEP VCF (SNV session only).
        germline_vcf: Germline VEP VCF (SNV session only).
        targets: Capture targets BED.
    """

    tumor_bam: Path | None = None
    normal_bam: Path | None = None
    tumor_evidence_bam: Path | None = None
    normal_evidence_bam: Path | None = None
    tumor_segments: Path | None = None
    tumor_bins: Path | None = None
    normal_segments: Path | None = None
    normal_bins: Path | None = None
    germline_afs: Path | None = None
    somatic_vcf: Path | None = None
    germline_vcf: Path | None = None
    targets: Path | None = None


@dataclass(frozen=True)
class TrackSpec:
    """One input role placed in a session.

    Attributes:
        role: ``SessionInputs`` attribute holding the file.
        kind: Track layout.
        name: Track name; ``None`` uses the file name.
    """

    role: str
    kind: TrackKind
    name: str | None = None


@dataclass(frozen=True)
class SessionLayout:
    """Template, track order and alignment panel size for one session file.

    Attributes:
        template: Bundled template file name.
        tracks: Tracks in display order.
        panel_height: Height of each alignment panel.
        panel_width: Width of each alignment panel.
    """

    template: str
    tracks: tuple[TrackSpec, ...]
    panel_height: int
    panel_width: int


_CNV_AND_AFS_TRACKS = (
    TrackSpec("normal_segments", TrackKind.SEGMENTS, "gDNA_segments"),
    TrackSpec("normal_bins", TrackKind.BINS, "gDNA_bins"),
    TrackSpec("tumor_segments", TrackKind.SEGMENTS, "tDNA_segments"),
    TrackSpec("tumor_bins", TrackKind.BINS, "tDNA_bins"),
    TrackSpec("germline_afs", TrackKind.AFS, "SNP_tDNA"),
)

# Track order follows the legacy script.
SNV_LAYOUT = SessionLayout(
    template=_SNV_TEMPLATE,
    tracks=(
        TrackSpec("tumor_bam", TrackKind.ALIGNMENT),
        TrackSpec("normal_bam", TrackKind.ALIGNMENT),
        *_CNV_AND_AFS_TRACKS,
        TrackSpec("germline_vcf", TrackKind.VCF, "GERMLINE VEP VCF"),
        TrackSpec("somatic_vcf", TrackKind.VCF, "SOMATIC VEP VCF"),
    ),
    panel_height=250,
    panel_width=1800,
)

SV_LAYOUT = SessionLayout(
    template=_SV_TEMPLATE,
    tracks=(
        TrackSpec("tumor_bam", TrackKind.ALIGNMENT),
        TrackSpec("normal_bam", TrackKind.ALIGNMENT),
        *_CNV_AND_AFS_TRACKS,
        TrackSpec("normal_evidence_bam", TrackKind.ALIGNMENT),
        TrackSpec("tumor_evidence_bam", TrackKind.ALIGNMENT),
    ),
    panel_height=50,
    panel_width=2543,
)


@dataclass
class _SessionParts:
    resources: list[str] = field(default_factory=list)
    data_tracks: list[str] = field(default_factory=list)
    alignment_panels: list[str] = field(default_factory=list)


# -- inputs -----------------------------------------------------------------


def drop_empty_inputs(inputs: SessionInputs) -> SessionInputs:
    """Return ``inputs`` with empty files set to ``None``, warning for each."""
    kept: dict[str, Path | None] = {}
    for f in fields(inputs):
        path: Path | None = getattr(inputs, f.name)
        if path is not None and path.stat().st_size == 0:
            logger.warning("Skipping empty file for --%s: %s", f.name.replace("_", "-"), path)
            path = None
        kept[f.name] = path
    return SessionInputs(**kept)


def _has_index(bam: Path) -> bool:
    return bam.with_name(bam.name + ".bai").exists() or bam.with_suffix(".bai").exists()


def _warn_missing_bam_indexes(inputs: SessionInputs) -> None:
    for bam in (
        inputs.tumor_bam,
        inputs.normal_bam,
        inputs.tumor_evidence_bam,
        inputs.normal_evidence_bam,
    ):
        if bam is not None and not _has_index(bam):
            logger.warning("No .bai index next to %s; IGV cannot load it without one", bam)


def _warn_missing_reference_dir(reference_dir: Path) -> None:
    if not reference_dir.is_dir():
        logger.warning("Reference directory %s does not exist on this host", reference_dir)


# -- rendering --------------------------------------------------------------


def _attr(value: str) -> str:
    """Escape a value for a double-quoted XML attribute."""
    return escape(value, {'"': "&quot;"})


def format_path(path: Path, outdir: Path, style: PathStyle) -> str:
    """Render an input path as absolute or relative to ``outdir``."""
    absolute = path.resolve()
    if style is PathStyle.RELATIVE:
        return os.path.relpath(absolute, outdir.resolve())
    return str(absolute)


def bam_color(bam_name: str) -> str:
    """Return the alignment track colour for a BAM file name."""
    match = _BAM_SV_TYPE.search(bam_name)
    return _SV_TYPE_COLORS[match.group(2)] if match else _DEFAULT_BAM_COLOR


def _add_track(
    parts: _SessionParts, spec: TrackSpec, file: Path, path: str, layout: SessionLayout
) -> None:
    parts.resources.append(_RESOURCE.format(path=_attr(path)))
    if spec.kind is TrackKind.ALIGNMENT:
        tracks = _ALIGNMENT_TRACKS.format(
            path=_attr(path), name=_attr(file.name), color=bam_color(file.name)
        )
        parts.alignment_panels.append(
            _PANEL.format(
                height=layout.panel_height,
                width=layout.panel_width,
                name=_attr(file.name),
                tracks=tracks,
            )
        )
        return
    name = spec.name or file.name
    parts.data_tracks.append(
        _DATA_TRACK_TEMPLATES[spec.kind].format(path=_attr(path), name=_attr(name))
    )


def render_session(
    template: str,
    layout: SessionLayout,
    inputs: SessionInputs,
    *,
    outdir: Path,
    path_style: PathStyle,
    genome: str,
    reference_dir: Path,
) -> str:
    """Fill a session template with the tracks of ``layout``.

    Args:
        template: Template text with ``{resources}``, ``{data_tracks}``,
            ``{alignment_panels}``, ``{capture_bed_track}``, ``{genome}`` and
            ``{reference_dir}`` placeholders.
        layout: Track order and alignment panel size.
        inputs: Input files; ``None`` roles are left out.
        outdir: Directory the session is written to (base for relative paths).
        path_style: Absolute or relative input paths.
        genome: IGV genome ID.
        reference_dir: Directory holding the template's reference tracks.

    Raises:
        InputValidationError: If the template has unknown or malformed placeholders.
    """
    parts = _SessionParts()
    for spec in layout.tracks:
        file: Path | None = getattr(inputs, spec.role)
        if file is not None:
            _add_track(parts, spec, file, format_path(file, outdir, path_style), layout)

    capture_bed_track = ""
    if inputs.targets is not None:
        path = _attr(format_path(inputs.targets, outdir, path_style))
        parts.resources.insert(0, _RESOURCE.format(path=path))
        capture_bed_track = _CAPTURE_BED_TRACK.format(path=path, name=_attr(inputs.targets.name))

    try:
        return template.format(
            resources="\n        ".join(parts.resources),
            data_tracks="\n        ".join(parts.data_tracks),
            alignment_panels="\n    ".join(parts.alignment_panels),
            capture_bed_track=capture_bed_track,
            genome=_attr(genome),
            reference_dir=_attr(str(reference_dir.resolve())),
        )
    except (KeyError, IndexError, ValueError) as exc:
        raise InputValidationError(f"Session template has an invalid placeholder: {exc}") from exc


def _read_template(custom: Path | None, bundled: str) -> str:
    if custom is None:
        return _TEMPLATES.joinpath(bundled).read_text(encoding="utf-8")
    try:
        return custom.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise InputValidationError(f"Could not read template {custom}: {exc}") from exc


# -- CLI --------------------------------------------------------------------


def _file_option(flag: str, help_text: str) -> Any:
    return typer.Option(flag, exists=True, dir_okay=False, help=help_text)


def igv_session(
    outdir: Annotated[Path, typer.Option("--outdir", "-o", help="Directory for the sessions.")],
    reference_dir: Annotated[
        Path,
        typer.Option(
            "--reference-dir", file_okay=False, help="Directory with the reference tracks."
        ),
    ],
    genome: Annotated[str, typer.Option("--genome", help="IGV genome ID.")] = "hg19",
    tumor_bam: Annotated[Path | None, _file_option("--tumor-bam", "Tumor BAM.")] = None,
    normal_bam: Annotated[Path | None, _file_option("--normal-bam", "Normal BAM.")] = None,
    tumor_evidence_bam: Annotated[
        Path | None,
        _file_option("--tumor-evidence-bam", "Tumor GRIDSS evidence BAM (SV session)."),
    ] = None,
    normal_evidence_bam: Annotated[
        Path | None,
        _file_option("--normal-evidence-bam", "Normal GRIDSS evidence BAM (SV session)."),
    ] = None,
    tumor_segments: Annotated[
        Path | None, _file_option("--tumor-segments", "Tumor CNV segments bedGraph.")
    ] = None,
    tumor_bins: Annotated[
        Path | None, _file_option("--tumor-bins", "Tumor CNV bins bedGraph.")
    ] = None,
    normal_segments: Annotated[
        Path | None, _file_option("--normal-segments", "Normal CNV segments bedGraph.")
    ] = None,
    normal_bins: Annotated[
        Path | None, _file_option("--normal-bins", "Normal CNV bins bedGraph.")
    ] = None,
    germline_afs: Annotated[
        Path | None,
        _file_option("--germline-afs", "Germline SNP allele frequencies bedGraph."),
    ] = None,
    somatic_vcf: Annotated[
        Path | None, _file_option("--somatic-vcf", "Somatic VEP VCF (SNV session).")
    ] = None,
    germline_vcf: Annotated[
        Path | None, _file_option("--germline-vcf", "Germline VEP VCF (SNV session).")
    ] = None,
    targets: Annotated[Path | None, _file_option("--targets", "Capture targets BED.")] = None,
    path_style: Annotated[
        PathStyle,
        typer.Option("--path-style", help="Write input paths absolute or relative to --outdir."),
    ] = PathStyle.ABSOLUTE,
    snv_template: Annotated[
        Path | None, _file_option("--snv-template", "Override the bundled SNV template.")
    ] = None,
    sv_template: Annotated[
        Path | None, _file_option("--sv-template", "Override the bundled SV template.")
    ] = None,
    dry_run: Annotated[
        bool, typer.Option("--dry-run", help="Validate inputs without writing sessions.")
    ] = False,
    verbose: VerboseOpt = False,
    quiet: QuietOpt = False,
) -> None:
    """Write IGV session files (SNV and SV) for one sample's results."""
    configure_logging(verbose=verbose, quiet=quiet)
    inputs = SessionInputs(
        tumor_bam=tumor_bam,
        normal_bam=normal_bam,
        tumor_evidence_bam=tumor_evidence_bam,
        normal_evidence_bam=normal_evidence_bam,
        tumor_segments=tumor_segments,
        tumor_bins=tumor_bins,
        normal_segments=normal_segments,
        normal_bins=normal_bins,
        germline_afs=germline_afs,
        somatic_vcf=somatic_vcf,
        germline_vcf=germline_vcf,
        targets=targets,
    )
    try:
        if all(getattr(inputs, f.name) is None for f in fields(inputs)):
            raise InputValidationError("No input files given.")
        inputs = drop_empty_inputs(inputs)
        _warn_missing_bam_indexes(inputs)
        _warn_missing_reference_dir(reference_dir)

        sessions = {
            SNV_SESSION_NAME: (SNV_LAYOUT, _read_template(snv_template, SNV_LAYOUT.template)),
            SV_SESSION_NAME: (SV_LAYOUT, _read_template(sv_template, SV_LAYOUT.template)),
        }
        rendered = {
            name: render_session(
                template,
                layout,
                inputs,
                outdir=outdir,
                path_style=path_style,
                genome=genome,
                reference_dir=reference_dir,
            )
            for name, (layout, template) in sessions.items()
        }
    except Vcf2CuratorError as exc:
        fail(exc)

    if dry_run:
        logger.info("Dry run: would write %s to %s", ", ".join(rendered), outdir)
        return

    outdir.mkdir(parents=True, exist_ok=True)
    for name, xml in rendered.items():
        (outdir / name).write_text(xml, encoding="utf-8")
        logger.info("Wrote %s", outdir / name)
