"""``svs`` command: GRIDSS structural variant VCF -> Curator table.

Replaces the legacy ``gridss_svannotate.R``, ``generate_evidence_bam.py`` and
the GRIDSS path of ``generateIGVnavInput_SV.py``. Steps:

1. Read PASS breakend records and classify each event (:mod:`vcf2curator.utils.gridss`).
2. Optionally (``--bam``) drop events failing the SAME_START_READS filter and
   write an evidence BAM of the remaining events' supporting reads.
3. Annotate genes, CGC, the CURATOR flag and curation tags, then write the table.

Behavioural notes vs. the legacy scripts (intentional, flagged for review):

* INS ``SV_LENGTH`` is the number of inserted bases; legacy counted every letter
  in ALT, including the reference base and any ``chr`` prefix.
* ``CGC_ANN`` values are sorted; legacy joined a set, so the order varied.
* ``GENE_A``/``GENE_B`` list each overlapping gene once.
* The exon check spans the two breakends whichever comes first.
* The evidence BAM holds only reads of events written to the table (legacy took
  every record), each alignment once.
* The TOOL column, svcaller/svaba inputs, ``.mut`` files and IGV colour types
  are dropped.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import pysam
import typer

from vcf2curator.annotations.cgc import load_cgc_by_symbol
from vcf2curator.annotations.curation import annotate_curate_tag, load_curation
from vcf2curator.annotations.intervals import (
    IntervalIndex,
    load_bed,
    near_any,
    overlapping_names,
    overlaps_any,
)
from vcf2curator.commands._shared import QuietOpt, VerboseOpt, fail
from vcf2curator.curator.schema import SV_COLUMNS, VcfType
from vcf2curator.curator.writer import CuratorWriter, Row
from vcf2curator.exceptions import InputValidationError, Vcf2CuratorError
from vcf2curator.logging_config import configure_logging, get_logger
from vcf2curator.utils.evidence import (
    AlignmentsByName,
    collect_alignments,
    has_same_start_reads,
    write_evidence_bam,
)
from vcf2curator.utils.gridss import Breakend, SvType, parse_breakend, simple_type

logger = get_logger()

_VALID_CHROMS = frozenset([str(i) for i in range(1, 23)] + ["X", "Y"])


@dataclass(frozen=True)
class SvFilterConfig:
    """Tunable filter thresholds, all defaulting to the legacy behaviour.

    Attributes:
        germline_min_support: Germline events with fewer supporting reads are
            set to CURATOR=NO.
        target_padding: Distance (bp) around a target that still counts as on target.
        same_start_fraction: SAME_START_READS fails when one read signature
            reaches this share of the supporting reads.
        min_distinct_reads: SAME_START_READS fails below this many distinct
            read signatures.
    """

    germline_min_support: int = 40
    target_padding: int = 150
    same_start_fraction: float = 0.6
    min_distinct_reads: int = 3


@dataclass(frozen=True)
class SvAnnotations:
    """Annotation resources used to build SV rows.

    Attributes:
        genes: Gene intervals (BED name = gene symbol).
        exons: Exon intervals.
        cgc: Gene symbol -> CGC annotation.
        targets: Capture target intervals, or ``None`` for whole-genome data.
        curation: Curation table, or ``None`` when not supplied.
    """

    genes: IntervalIndex
    exons: IntervalIndex
    cgc: dict[str, str]
    targets: IntervalIndex | None
    curation: pd.DataFrame | None


@dataclass(frozen=True)
class SvEvent:
    """One breakend record kept for output.

    Attributes:
        bnd: The parsed breakend and its partner.
        svtype: Simple event type.
        support: Supporting read count (INFO/VF).
        vaf: Variant allele frequency as written to the table.
        read_names: Supporting read names (INFO/BPNAMES), empty if absent.
    """

    bnd: Breakend
    svtype: SvType
    support: int
    vaf: str
    read_names: tuple[str, ...]


# -- reading ----------------------------------------------------------------


def _is_valid_chrom(chrom: str) -> bool:
    return chrom.removeprefix("chr") in _VALID_CHROMS


def _is_pass(record: pysam.VariantRecord) -> bool:
    keys = list(record.filter.keys())
    return not keys or keys == ["PASS"]


def _format_vaf(value: Any) -> str:
    """Render a VAF value (scalar or tuple), trimming float32 widening noise."""
    values = value if isinstance(value, tuple) else (value,)
    return ",".join(str(round(v, 4)) if isinstance(v, float) else str(v) for v in values)


def _check_header(vcf: pysam.VariantFile, vcftype: VcfType, need_read_names: bool) -> None:
    required = ["EVENT", "VF"]
    if vcftype is VcfType.SOMATIC:
        required.append("TAF")
    if need_read_names:
        required.append("BPNAMES")
    missing = [key for key in required if key not in vcf.header.info]
    if missing:
        raise InputValidationError(
            f"GRIDSS VCF header is missing INFO fields {missing}; is this a GRIDSS VCF?"
        )
    if vcftype is VcfType.GERMLINE and "AF" not in vcf.header.formats:
        raise InputValidationError("Germline GRIDSS VCF header is missing FORMAT/AF.")


def _vaf(record: pysam.VariantRecord, vcftype: VcfType) -> str:
    if vcftype is VcfType.SOMATIC:
        return _format_vaf(record.info["TAF"])
    return _format_vaf(record.samples[0]["AF"])


def read_events(vcf_path: str | Path, vcftype: VcfType, need_read_names: bool) -> list[SvEvent]:
    """Read the PASS events to report from a GRIDSS VCF.

    Keeps one record per EVENT, except translocations, where both breakend
    records are kept (one row per side). Single breakends and events touching a
    chromosome outside 1-22/X/Y (with or without ``chr``) are skipped.

    Args:
        vcf_path: GRIDSS VCF path.
        vcftype: Somatic or germline call set.
        need_read_names: Require INFO/BPNAMES (needed for the BAM steps).

    Raises:
        InputValidationError: If the VCF cannot be read or lacks GRIDSS fields.
    """
    try:
        vcf = pysam.VariantFile(str(vcf_path))
    except (OSError, ValueError) as exc:
        raise InputValidationError(f"Could not open VCF {vcf_path}: {exc}") from exc

    events: list[SvEvent] = []
    seen_events: set[str] = set()
    with vcf:
        _check_header(vcf, vcftype, need_read_names)
        for record in vcf:
            if not _is_pass(record) or not record.alts or record.ref is None:
                continue
            bnd = parse_breakend(record.chrom, record.pos, record.ref, record.alts[0])
            if bnd is None:
                continue
            if not (_is_valid_chrom(bnd.chrom) and _is_valid_chrom(bnd.partner_chrom)):
                continue
            event_id = record.info.get("EVENT")
            if event_id is None:
                raise InputValidationError(
                    f"{record.chrom}:{record.pos} ({record.id}) has no INFO/EVENT"
                )
            svtype = simple_type(bnd)
            if event_id in seen_events:
                continue
            if svtype is not SvType.TRA:
                seen_events.add(event_id)
            read_names = record.info.get("BPNAMES", ())
            if need_read_names and not read_names:
                raise InputValidationError(
                    f"{record.chrom}:{record.pos} ({record.id}) has no INFO/BPNAMES"
                )
            events.append(
                SvEvent(
                    bnd=bnd,
                    svtype=svtype,
                    support=int(record.info["VF"]),
                    vaf=_vaf(record, vcftype),
                    read_names=tuple(read_names),
                )
            )
    logger.info("Read %d SV events from %s", len(events), vcf_path)
    return events


# -- annotation -------------------------------------------------------------


def _genes(annotations: SvAnnotations, chrom: str, start: int, end: int) -> str:
    names = overlapping_names(annotations.genes, chrom, start, end)
    return ",".join(names) if names else "NA"


def _cgc_annotations(cgc: dict[str, str], *gene_fields: str) -> list[str]:
    genes = {g for field in gene_fields for g in field.split(",")}
    return sorted({cgc[g] for g in genes if cgc.get(g)})


def _curate(annotations: SvAnnotations, genes_sorted: str, project: str, vcftype: VcfType) -> str:
    if annotations.curation is None:
        return ""
    tags = [
        annotate_curate_tag(annotations.curation, project, gene, vcftype.value)
        for gene in genes_sorted.split(",")
    ]
    return ",".join(tag for tag in tags if tag)


def _curator_flag(
    event: SvEvent,
    coords: dict[str, Any],
    has_cgc: bool,
    annotations: SvAnnotations,
    vcftype: VcfType,
    filters: SvFilterConfig,
) -> str:
    """Decide CURATOR: base rule, then germline support, then exon overlap."""
    chrom_a, chrom_b = coords["CHROM_A"], coords["CHROM_B"]
    start_a, end_a = coords["START_A"], coords["END_A"]
    start_b, end_b = coords["START_B"], coords["END_B"]
    pad = filters.target_padding

    if annotations.targets is None:
        on = has_cgc
    elif event.svtype is SvType.TRA:
        on = near_any(annotations.targets, chrom_a, (start_a, end_a), pad) or near_any(
            annotations.targets, chrom_b, (start_b, end_b), pad
        )
    else:
        on = near_any(annotations.targets, chrom_a, (start_a, end_b), pad)

    if vcftype is VcfType.GERMLINE and event.support < filters.germline_min_support:
        on = False
    if on and event.svtype is not SvType.TRA:
        span_start, span_end = min(start_a, end_b), max(start_a, end_b)
        on = overlaps_any(annotations.exons, chrom_a, span_start, span_end)
    return "YES" if on else "NO"


def build_row(
    event: SvEvent,
    annotations: SvAnnotations,
    *,
    vcftype: VcfType,
    sdid: str,
    project: str,
    filters: SvFilterConfig,
) -> Row:
    """Build the Curator row for one event, using the legacy coordinate scheme."""
    bnd = event.bnd
    if event.svtype is SvType.TRA:
        start_a, end_a = bnd.pos, bnd.pos + 1
    else:
        start_a, end_a = bnd.pos - 1, bnd.pos
    start_b, end_b = bnd.partner_pos - 1, bnd.partner_pos
    coords: dict[str, Any] = {
        "CHROM_A": bnd.chrom,
        "START_A": start_a,
        "END_A": end_a,
        "CHROM_B": bnd.partner_chrom,
        "START_B": start_b,
        "END_B": end_b,
    }

    if event.svtype is SvType.INS:
        length: Any = len(bnd.inserted)
    elif bnd.chrom == bnd.partner_chrom:
        length = abs(end_b - start_a)
    else:
        length = "NA"

    gene_a = _genes(annotations, bnd.chrom, start_a, end_a)
    gene_b = _genes(annotations, bnd.partner_chrom, start_b, end_b)
    genes_sorted = ",".join(sorted([gene_a, gene_b]))
    cgc = _cgc_annotations(annotations.cgc, gene_a, gene_b)

    return {
        **coords,
        "IGV_COORD": f"{bnd.chrom}:{bnd.pos} {bnd.partner_chrom}:{start_b}",
        "SVTYPE": event.svtype.value,
        "SV_LENGTH": length,
        "SUPPORT_READS": event.support,
        "VAF": event.vaf,
        "SDID": sdid,
        "SAMPLE": vcftype.value,
        "GENE_A": gene_a,
        "GENE_B": gene_b,
        "GENE_A-GENE_B-sorted": genes_sorted,
        "CGC_ANN": ",".join(cgc),
        "CURATOR": _curator_flag(event, coords, bool(cgc), annotations, vcftype, filters),
        "CURATE": _curate(annotations, genes_sorted, project, vcftype),
    }


def _sort_key(row: Row) -> tuple[Any, ...]:
    return (
        row["GENE_A-GENE_B-sorted"],
        row["CHROM_A"],
        row["START_A"],
        row["CHROM_B"],
        row["START_B"],
    )


def drop_same_start_events(
    events: list[SvEvent], alignments: AlignmentsByName, filters: SvFilterConfig
) -> list[SvEvent]:
    """Remove events whose supporting reads fail the SAME_START_READS filter."""
    kept = [
        e
        for e in events
        if not has_same_start_reads(
            e.read_names,
            alignments,
            fraction=filters.same_start_fraction,
            min_distinct=filters.min_distinct_reads,
        )
    ]
    logger.info("SAME_START_READS removed %d of %d events", len(events) - len(kept), len(events))
    return kept


def build_rows(
    events: list[SvEvent],
    annotations: SvAnnotations,
    *,
    vcftype: VcfType,
    sdid: str,
    project: str,
    filters: SvFilterConfig,
) -> Iterator[Row]:
    """Yield Curator rows for ``events`` in legacy sort order."""
    rows = [
        build_row(e, annotations, vcftype=vcftype, sdid=sdid, project=project, filters=filters)
        for e in events
    ]
    yield from sorted(rows, key=_sort_key)


# -- CLI command -------------------------------------------------------------


def _check_options(
    targets: Path | None, wgs: bool, bam: Path | None, evidence: Path | None
) -> None:
    if (targets is None) == (not wgs):
        raise InputValidationError("Give exactly one of --targets (panel) or --wgs.")
    if evidence is not None and bam is None:
        raise InputValidationError("--evidence-bam requires --bam.")


def svs(
    input_file: Annotated[
        Path,
        typer.Option("--input", "-i", exists=True, dir_okay=False, help="GRIDSS SV VCF."),
    ],
    output: Annotated[Path, typer.Option("--output", "-o", help="Destination Curator TSV.")],
    vcftype: Annotated[VcfType, typer.Option("--vcftype", help="Call set type.")],
    sdid: Annotated[str, typer.Option("--sdid", help="Sample ID written to the SDID column.")],
    project: Annotated[
        str, typer.Option("--project", help="Project prefix used for curation lookups.")
    ],
    genes: Annotated[
        Path,
        typer.Option(
            "--genes", exists=True, dir_okay=False, help="Genes BED (chrom, start, end, symbol)."
        ),
    ],
    exons: Annotated[Path, typer.Option("--exons", exists=True, dir_okay=False, help="Exons BED.")],
    cgc: Annotated[
        Path, typer.Option("--cgc", exists=True, dir_okay=False, help="Cancer Gene Census TSV.")
    ],
    targets: Annotated[
        Path | None,
        typer.Option(
            "--targets", exists=True, dir_okay=False, help="Capture targets BED (panel data)."
        ),
    ] = None,
    wgs: Annotated[
        bool, typer.Option("--wgs", help="Whole-genome data: CURATOR uses CGC genes.")
    ] = False,
    curation: Annotated[
        Path | None,
        typer.Option(
            "--curation-ann", exists=True, dir_okay=False, help="Curation annotation CSV."
        ),
    ] = None,
    bam: Annotated[
        Path | None,
        typer.Option(
            "--bam",
            exists=True,
            dir_okay=False,
            help="BAM the SVs were called from; enables the SAME_START_READS filter.",
        ),
    ] = None,
    evidence_bam: Annotated[
        Path | None,
        typer.Option("--evidence-bam", help="Write supporting reads to this BAM (needs --bam)."),
    ] = None,
    same_start_fraction: Annotated[
        float,
        typer.Option("--same-start-fraction", help="SAME_START_READS: max share of one signature."),
    ] = 0.6,
    min_distinct_reads: Annotated[
        int,
        typer.Option("--min-distinct-reads", help="SAME_START_READS: min distinct signatures."),
    ] = 3,
    germline_min_support: Annotated[
        int,
        typer.Option("--germline-min-support", help="Germline: min supporting reads for CURATOR."),
    ] = 40,
    target_padding: Annotated[
        int, typer.Option("--target-padding", help="Bases around targets counted as on target.")
    ] = 150,
    dry_run: Annotated[
        bool, typer.Option("--dry-run", help="Parse and report counts without writing output.")
    ] = False,
    verbose: VerboseOpt = False,
    quiet: QuietOpt = False,
) -> None:
    """Convert a GRIDSS structural variant VCF into a Curator table."""
    logger = configure_logging(verbose=verbose, quiet=quiet)
    filters = SvFilterConfig(
        germline_min_support=germline_min_support,
        target_padding=target_padding,
        same_start_fraction=same_start_fraction,
        min_distinct_reads=min_distinct_reads,
    )

    try:
        _check_options(targets, wgs, bam, evidence_bam)
        annotations = SvAnnotations(
            genes=load_bed(genes, require_name=True),
            exons=load_bed(exons),
            cgc=load_cgc_by_symbol(cgc),
            targets=load_bed(targets) if targets else None,
            curation=load_curation(curation) if curation else None,
        )
        events = read_events(input_file, vcftype, need_read_names=bam is not None)
        alignments: AlignmentsByName = {}
        if bam is not None:
            wanted = {name for e in events for name in e.read_names}
            alignments = collect_alignments(bam, wanted)
            events = drop_same_start_events(events, alignments, filters)
        rows = build_rows(
            events, annotations, vcftype=vcftype, sdid=sdid, project=project, filters=filters
        )

        if dry_run:
            count = sum(1 for _ in rows)
            logger.info("Dry run: %d rows would be written to %s", count, output)
            typer.echo(f"[dry-run] {count} rows")
            return
        count = CuratorWriter(SV_COLUMNS).write(rows, output)
        typer.echo(f"Wrote {count} rows to {output}")
        if bam is not None and evidence_bam is not None:
            names = {name for e in events for name in e.read_names}
            write_evidence_bam(bam, alignments, names, evidence_bam)
            typer.echo(f"Wrote evidence BAM {evidence_bam}")
    except Vcf2CuratorError as exc:
        fail(exc)
