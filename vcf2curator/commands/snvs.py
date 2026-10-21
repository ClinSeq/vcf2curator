"""``snvs`` command: SNV/indel VCF -> Curator table.

Port of the legacy ``generateIGVnavInput.py``. Reads a VEP-annotated somatic or
germline VCF and emits the IGVNav/Curator table. The VCF backend is ``pysam``
(the legacy script used the abandoned PyVCF); the transformation, filters and
output columns are preserved.

Behavioural notes vs. the original script (intentional, flagged for review):

* The ``ALT`` column contains the alt allele(s) joined by ``,`` (e.g. ``T``)
  rather than PyVCF's Python list repr (e.g. ``[T]``).
* ``T_VAF`` / ``N_VAF`` from FORMAT ``AF``/``VAF`` are scalarised (first value)
  and rounded to 4 dp to remove pysam's float32 widening noise, instead of
  being written as a one-element list.
* Somatic tumor/normal samples are identified by the ``-N-`` naming convention;
  if neither can be found the parser fails loudly instead of silently using
  broken placeholder sample names.
* Hotspot and CGC annotations are skipped (left blank) when their resource files
  are not supplied, rather than raising; supply the files to enable them.
* The IGVnav symlink side-effect and the unused ``--vardict`` input are dropped.
"""

from __future__ import annotations

import os
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Annotated, Any

import pysam
import typer

from vcf2curator.annotations.aminoacid import tri_to_single_aa
from vcf2curator.annotations.cgc import CgcLookup, load_cgc
from vcf2curator.annotations.curation import annotate_curate_tag, load_curation
from vcf2curator.annotations.hotspot import (
    annotate_hotspot,
    check_noncoding_hotspot,
    load_hotspot,
)
from vcf2curator.annotations.oncokb import OncoKBLookup, load_oncokb, lookup_oncogenicity
from vcf2curator.commands._shared import QuietOpt, VerboseOpt, fail
from vcf2curator.curator.schema import VcfType, columns_for
from vcf2curator.curator.writer import CuratorWriter, Row
from vcf2curator.exceptions import InputValidationError, Vcf2CuratorError
from vcf2curator.logging_config import configure_logging, get_logger

logger = get_logger()

# Default FILTER values kept, per VCF type (legacy behaviour).
_DEFAULT_SOMATIC_FILTERS = frozenset({"PASS", "LowQual"})
_DEFAULT_GERMLINE_FILTERS = frozenset({"PASS"})


@dataclass(frozen=True)
class SnvFilterConfig:
    """Tunable filter thresholds, all defaulting to the legacy behaviour.

    Attributes:
        impacts: VEP IMPACT values considered actionable.
        keep_filters: FILTER values to keep; ``None`` selects the per-VCF-type
            default (somatic: PASS/LowQual, germline: PASS).
        include_splice: Keep splice-region variants regardless of IMPACT.
        include_hotspot: Keep hotspot variants regardless of IMPACT (somatic).
        wgs_min_alt: Minimum tumor ALT read support for WGS somatic calls.
        exclude_hla: Drop germline WGS variants in HLA genes.
        require_pathogenic_missense: Drop germline missense variants unless the
            clinical significance contains "pathogenic".
        require_pathogenic_splice: Drop germline splice variants unless the
            clinical significance contains "pathogenic".
    """

    impacts: frozenset[str] = frozenset({"HIGH", "MODERATE"})
    keep_filters: frozenset[str] | None = None
    include_splice: bool = True
    include_hotspot: bool = True
    wgs_min_alt: int = 5
    exclude_hla: bool = True
    require_pathogenic_missense: bool = True
    require_pathogenic_splice: bool = True

    def resolved_filters(self, vcftype: VcfType) -> frozenset[str]:
        """Return the FILTER allow-list, falling back to the per-type default."""
        if self.keep_filters is not None:
            return self.keep_filters
        return _DEFAULT_SOMATIC_FILTERS if vcftype is VcfType.SOMATIC else _DEFAULT_GERMLINE_FILTERS


@dataclass
class SnvParser:
    """Parse a VEP-annotated SNV/indel VCF into Curator rows.

    Attributes:
        vcftype: Somatic or germline call set.
        oncokb_path: Path to the OncoKB all-variants TSV (required).
        cgc_path: Optional Cancer Gene Census TSV.
        hotspot_snv_path: Optional cancer hotspot SNV CSV.
        hotspot_indel_path: Optional cancer hotspot indel CSV.
        curation_path: Optional curation annotation CSV.
        wgs: Apply whole-genome-sequencing filters instead of panel filters.
        filters: Tunable filter thresholds.
    """

    vcftype: VcfType
    oncokb_path: str | Path
    cgc_path: str | Path | None = None
    hotspot_snv_path: str | Path | None = None
    hotspot_indel_path: str | Path | None = None
    curation_path: str | Path | None = None
    wgs: bool = False
    filters: SnvFilterConfig = field(default_factory=SnvFilterConfig)

    def columns(self) -> tuple[str, ...]:
        """Return the Curator columns for this VCF type."""
        return columns_for(self.vcftype)

    # -- resource loading ---------------------------------------------------

    def _load_resources(self) -> _Resources:
        oncokb = load_oncokb(self.oncokb_path)
        cgc: CgcLookup = load_cgc(self.cgc_path) if self.cgc_path else {}
        hotspot_snv = load_hotspot(self.hotspot_snv_path, "snv") if self.hotspot_snv_path else None
        hotspot_indel = (
            load_hotspot(self.hotspot_indel_path, "indel") if self.hotspot_indel_path else None
        )
        curation = load_curation(self.curation_path) if self.curation_path else None
        return _Resources(oncokb, cgc, hotspot_snv, hotspot_indel, curation)

    # -- public API ---------------------------------------------------------

    def _open(self, input_path: str | Path) -> pysam.VariantFile:
        try:
            return pysam.VariantFile(str(input_path))
        except (OSError, ValueError) as exc:
            raise InputValidationError(f"Could not open VCF {input_path}: {exc}") from exc

    def _csq_keys(self, vcf: pysam.VariantFile) -> list[str]:
        if "CSQ" not in vcf.header.info:
            raise InputValidationError(
                "Input VCF has no CSQ INFO field; a VEP-annotated VCF is required."
            )
        description = vcf.header.info["CSQ"].description
        if not description:
            raise InputValidationError("CSQ INFO header has no field description.")
        if "Format: " in description:
            fields = description.split("Format: ")[-1]
        else:
            fields = description.split(": ")[1]
        return fields.split("|")

    def validate(self, input_path: str | Path) -> None:
        """Check the VCF is readable, VEP-annotated, and has required CSQ keys."""
        vcf = self._open(input_path)
        try:
            csq_keys = self._csq_keys(vcf)
            required = {
                "SYMBOL",
                "Gene",
                "Amino_acids",
                "Protein_position",
                "IMPACT",
                "Consequence",
                "Feature",
                "HGVSc",
                "HGVSp",
                "CLIN_SIG",
                "BrcaEx_ClinicalSignificance",
            }
            missing = required - set(csq_keys)
            if missing:
                raise InputValidationError(
                    f"CSQ annotation is missing required fields: {sorted(missing)}"
                )
        finally:
            vcf.close()

    def parse(self, input_path: str | Path) -> Iterator[Row]:
        """Yield Curator rows parsed from a VEP-annotated VCF."""
        resources = self._load_resources()
        vcf = self._open(input_path)
        csq_keys = self._csq_keys(vcf)
        project = os.path.basename(str(input_path)).split("-")[0]
        keep_filters = self.filters.resolved_filters(self.vcftype)

        try:
            for record in vcf:
                if "CSQ" not in record.info:
                    continue
                row = self._process_record(record, csq_keys, resources, project, keep_filters)
                if row is not None:
                    yield row
        finally:
            vcf.close()

    # -- per-record processing ---------------------------------------------

    def _process_record(
        self,
        record: pysam.VariantRecord,
        csq_keys: list[str],
        resources: _Resources,
        project: str,
        keep_filters: frozenset[str],
    ) -> Row | None:
        can = _canonical_csq(record.info["CSQ"], csq_keys)
        try:
            gene = can["SYMBOL"]
            ensembl_id = can["Gene"]
            amino_acids = can["Amino_acids"].split("/")
            protein_positions = can["Protein_position"].split("/")
            clinsig = can["CLIN_SIG"]
            impact = can["IMPACT"]
            brca_ex = can["BrcaEx_ClinicalSignificance"]
            hgvsp = can["HGVSp"]
            consequence = can["Consequence"]
        except KeyError as exc:
            raise InputValidationError(f"CSQ annotation missing expected field {exc}") from exc

        gnomad = _gnomad_af(can, self.wgs)
        is_splice = "splice_region_variant" in consequence
        filter_col = _filter_col(record)
        alts = record.alts or ()
        ref = record.ref or ""
        is_snp = len(ref) == 1 and bool(alts) and all(len(a) == 1 for a in alts)

        oncogenicity = lookup_oncogenicity(resources.oncokb, gene, amino_acids, protein_positions)
        cgc_ann = resources.cgc.get(ensembl_id, ("", ""))[0] if resources.cgc else ""
        hotspot = self._annotate_hotspot(record, gene, hgvsp, consequence, is_snp, resources)
        curate = ""
        if resources.curation is not None:
            curate = annotate_curate_tag(resources.curation, project, gene, self.vcftype.value)
        aa_change = tri_to_single_aa(hgvsp.split("p.")[1]) if hgvsp else ""

        shared = _SharedFields(
            chrom=record.chrom,
            start=record.pos - 1,
            end=record.pos,
            ref=ref,
            alt=",".join(alts),
            gene=gene,
            ensembl_id=ensembl_id,
            impact=impact,
            consequence=consequence,
            transcript=can["Feature"],
            hgvsc=can["HGVSc"],
            hgvsp=hgvsp,
            amino_acid=aa_change,
            hotspot=hotspot,
            clinsig=clinsig,
            gnomad=gnomad,
            brca_ex=brca_ex,
            oncokb=oncogenicity,
            cgc_ann=cgc_ann,
            curate=curate,
        )

        if self.vcftype is VcfType.SOMATIC:
            return self._somatic_row(
                record,
                shared,
                impact,
                is_splice,
                hotspot,
                filter_col,
                keep_filters,
                can,
                alts,
            )
        return self._germline_row(
            record,
            shared,
            impact,
            is_splice,
            consequence,
            clinsig,
            filter_col,
            keep_filters,
            gene,
            cgc_ann,
            alts,
        )

    def _annotate_hotspot(
        self,
        record: pysam.VariantRecord,
        gene: str,
        hgvsp: str,
        consequence: str,
        is_snp: bool,
        resources: _Resources,
    ) -> str:
        hotspot = "Hotspot" if check_noncoding_hotspot(record.chrom, record.pos) else ""
        have_hotspot_tables = (
            resources.hotspot_snv is not None and resources.hotspot_indel is not None
        )
        if hgvsp and have_hotspot_tables:
            if is_snp:
                gpos = f"{record.chrom}:{record.pos}"
                hotspot = annotate_hotspot(
                    gene, hgvsp, "snv", consequence, resources.hotspot_snv, gpos
                )
            else:
                hotspot = annotate_hotspot(
                    gene, hgvsp, "indel", consequence, resources.hotspot_indel
                )
        return hotspot

    def _somatic_row(
        self,
        record: pysam.VariantRecord,
        shared: _SharedFields,
        impact: str,
        is_splice: bool,
        hotspot: str,
        filter_col: str,
        keep_filters: frozenset[str],
        can: dict[str, str],
        alts: tuple[str, ...],
    ) -> Row | None:
        normal, tumor = _somatic_samples(record)
        normal_dp, normal_alt, normal_vaf = _depth_alt_vaf(normal)
        tumor_dp, tumor_alt, tumor_vaf = _depth_alt_vaf(tumor)
        num_tools = _num_tools(record)
        rsid = can.get("Existing_variation", "")

        actionable = (
            impact in self.filters.impacts
            or (self.filters.include_splice and is_splice)
            or (self.filters.include_hotspot and hotspot == "Hotspot")
        )
        if not self.wgs:
            keep = filter_col in keep_filters and actionable
        else:
            keep = (
                tumor_alt is not None
                and tumor_alt >= self.filters.wgs_min_alt
                and (impact in self.filters.impacts)
            )
        if not keep:
            return None

        return {
            **shared.base(),
            "T_DP": tumor_dp,
            "T_ALT": tumor_alt,
            "T_VAF": tumor_vaf,
            "N_DP": normal_dp,
            "N_ALT": normal_alt,
            "N_VAF": normal_vaf,
            "RSID": rsid,
            "NUM_TOOLS": num_tools,
        }

    def _germline_row(
        self,
        record: pysam.VariantRecord,
        shared: _SharedFields,
        impact: str,
        is_splice: bool,
        consequence: str,
        clinsig: str,
        filter_col: str,
        keep_filters: frozenset[str],
        gene: str,
        cgc_ann: str,
        alts: tuple[str, ...],
    ) -> Row | None:
        names = list(record.samples)
        if not names:
            raise InputValidationError("Germline VCF record has no samples")
        normal = record.samples[names[0]]

        tumor_vaf: float = 0.0
        if len(names) > 1 and len(alts) == 1:
            tumor = record.samples[names[1]]
            t_ad = _get(tumor, "AD")
            t_dp = _get(tumor, "DP")
            if t_ad and t_dp:
                tumor_vaf = float(t_ad[1]) / float(t_dp)

        actionable = impact in self.filters.impacts or (self.filters.include_splice and is_splice)
        panel_ok = len(alts) == 1 and filter_col in keep_filters and actionable and not self.wgs
        wgs_ok = self.wgs and (
            (
                len(alts) == 1
                and impact == "HIGH"
                and not (self.filters.exclude_hla and gene.startswith("HLA"))
            )
            or ("pathogenic" in clinsig)
        )

        if panel_ok:
            if (
                "missense_variant" in consequence
                and self.filters.require_pathogenic_missense
                and "pathogenic" not in clinsig
            ):
                return None
            if is_splice and self.filters.require_pathogenic_splice and "pathogenic" not in clinsig:
                return None
        elif wgs_ok:
            if cgc_ann == "":
                return None
        else:
            return None

        n_dp = _get(normal, "DP")
        n_ad = _get(normal, "AD")
        if not (n_dp and n_ad):
            return None
        normal_alt = n_ad[1]
        normal_vaf = round(float(normal_alt) / float(n_dp), 2)

        return {
            **shared.base(),
            "N_DP": n_dp,
            "N_ALT": normal_alt,
            "N_VAF": normal_vaf,
            "T_VAF": round(tumor_vaf, 2),
            "RSID": record.id,
        }


# --------------------------------------------------------------------------
# Internal helpers and value objects
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class _Resources:
    """Bundle of loaded annotation resources for one parse run."""

    oncokb: OncoKBLookup
    cgc: CgcLookup
    hotspot_snv: Any
    hotspot_indel: Any
    curation: Any


@dataclass(frozen=True)
class _SharedFields:
    """Columns common to somatic and germline rows."""

    chrom: str
    start: int
    end: int
    ref: str
    alt: str
    gene: str
    ensembl_id: str
    impact: str
    consequence: str
    transcript: str
    hgvsc: str
    hgvsp: str
    amino_acid: str
    hotspot: str
    clinsig: str
    gnomad: str
    brca_ex: str
    oncokb: Any
    cgc_ann: str
    curate: str

    def base(self) -> dict[str, Any]:
        """Return the shared columns (CALL/TAG/NOTES intentionally blank)."""
        return {
            "CHROM": self.chrom,
            "START": self.start,
            "END": self.end,
            "REF": self.ref,
            "ALT": self.alt,
            "CALL": "",
            "TAG": "",
            "NOTES": "",
            "GENE": self.gene,
            "ENSEMBLID": self.ensembl_id,
            "IMPACT": self.impact,
            "CONSEQUENCE": self.consequence,
            "TRANSCRIPT": self.transcript,
            "HGVSc": self.hgvsc,
            "HGVSp": self.hgvsp,
            "AMINO_ACID": self.amino_acid,
            "HOTSPOT": self.hotspot,
            "CLIN_SIG": self.clinsig,
            "gnomAD": self.gnomad,
            "BRCAEx": self.brca_ex,
            "OncoKB": self.oncokb,
            "CGC_ANN": self.cgc_ann,
            "CURATE": self.curate,
        }


def _canonical_csq(csq_value: Any, csq_keys: list[str]) -> dict[str, str]:
    """Return the canonical transcript CSQ, or the first if none is canonical."""
    transcripts = (csq_value,) if isinstance(csq_value, str) else tuple(csq_value)
    parsed = [{csq_keys[idx]: ann for idx, ann in enumerate(t.split("|"))} for t in transcripts]
    for entry in parsed:
        if entry.get("CANONICAL") == "YES":
            return entry
    return parsed[0]


def _gnomad_af(can: dict[str, str], wgs: bool) -> str:
    """Select the gnomAD AF field matching the legacy precedence."""
    if "gnomAD_AF" in can:
        return can["gnomAD_AF"]
    return can["gnomADg_AF"] if wgs else can["gnomADe_AF"]


def _filter_col(record: pysam.VariantRecord) -> str:
    """Collapse the FILTER field to a single legacy-style value."""
    filters = [str(key) for key in record.filter.keys()]
    if not filters or filters == ["PASS"]:
        return "PASS"
    return filters[0]


def _get(sample: pysam.VariantRecordSample, key: str) -> Any:
    """Return a FORMAT value for a sample, or None if absent."""
    try:
        return sample[key]
    except KeyError:
        return None


def _has(sample: pysam.VariantRecordSample, key: str) -> bool:
    return key in sample.keys()


def _scalar(value: Any) -> Any:
    """Reduce a one-element tuple/list (e.g. Number=A field) to a scalar.

    Float values are rounded to 4 dp to strip the float32 -> float64 widening
    noise pysam introduces (e.g. 0.33 read back as 0.33000001311302185).
    """
    if isinstance(value, (tuple, list)):
        value = value[0] if value else None
    if isinstance(value, float):
        return round(value, 4)
    return value


def _depth_alt_vaf(sample: pysam.VariantRecordSample) -> tuple[Any, Any, Any]:
    """Extract (DP, ALT depth, VAF) with the legacy DP4 fallback."""
    if _has(sample, "DP"):
        dp = sample["DP"]
    else:
        dp = sum(sample["DP4"])
    if _has(sample, "AD"):
        alt = sample["AD"][1]
    else:
        alt = sum(sample["DP4"][2:])
    if _has(sample, "AF"):
        vaf = _scalar(sample["AF"])
    else:
        vaf = _scalar(_get(sample, "VAF"))
    return dp, alt, vaf


def _somatic_samples(
    record: pysam.VariantRecord,
) -> tuple[pysam.VariantRecordSample, pysam.VariantRecordSample]:
    """Identify the normal (``-N-``) and tumor samples by name convention."""
    names = [str(name) for name in record.samples]
    normal_name = next((n for n in names if "-N-" in n), None)
    tumor_name = next((n for n in names if "-N-" not in n), None)
    if normal_name is None or tumor_name is None:
        raise InputValidationError(
            "Could not identify tumor/normal samples by the '-N-' naming "
            f"convention in samples: {names}"
        )
    return record.samples[normal_name], record.samples[tumor_name]


def _num_tools(record: pysam.VariantRecord) -> int:
    """Return the NUM_TOOLS INFO value, defaulting to 1 when absent."""
    if "NUM_TOOLS" not in record.info:
        return 1
    value = record.info["NUM_TOOLS"]
    if isinstance(value, (tuple, list)):
        value = value[0]
    return int(value)


# -- CLI command -----------------------------------------------------------


def snvs(
    input_file: Annotated[
        Path,
        typer.Option(
            "--input",
            "-i",
            exists=True,
            dir_okay=False,
            help="VEP-annotated VCF (optionally bgzipped).",
        ),
    ],
    output: Annotated[Path, typer.Option("--output", "-o", help="Destination Curator TSV.")],
    vcftype: Annotated[VcfType, typer.Option("--vcftype", help="Call set type.")],
    oncokb: Annotated[
        Path, typer.Option("--oncokb", exists=True, dir_okay=False, help="OncoKB all-variants TSV.")
    ],
    cgc: Annotated[
        Path | None,
        typer.Option("--cgc", exists=True, dir_okay=False, help="Cancer Gene Census TSV."),
    ] = None,
    hotspot_snv: Annotated[
        Path | None,
        typer.Option("--hotspot-snv", exists=True, dir_okay=False, help="Cancer hotspot SNV CSV."),
    ] = None,
    hotspot_indel: Annotated[
        Path | None,
        typer.Option(
            "--hotspot-indel", exists=True, dir_okay=False, help="Cancer hotspot indel CSV."
        ),
    ] = None,
    curation: Annotated[
        Path | None,
        typer.Option(
            "--curation-ann", exists=True, dir_okay=False, help="Curation annotation CSV."
        ),
    ] = None,
    wgs: Annotated[bool, typer.Option("--wgs", help="Apply WGS filters.")] = False,
    impact: Annotated[
        list[str] | None,
        typer.Option(
            "--impact", help="Actionable IMPACT values (repeatable; default HIGH, MODERATE)."
        ),
    ] = None,
    keep_filter: Annotated[
        list[str] | None,
        typer.Option(
            "--keep-filter",
            help="FILTER values to keep (repeatable; default PASS/LowQual somatic, PASS germline).",
        ),
    ] = None,
    wgs_min_alt: Annotated[
        int, typer.Option("--wgs-min-alt", help="Min tumor ALT reads (WGS).")
    ] = 5,
    include_splice: Annotated[
        bool,
        typer.Option("--include-splice/--no-include-splice", help="Keep splice-region variants."),
    ] = True,
    include_hotspot: Annotated[
        bool, typer.Option("--include-hotspot/--no-include-hotspot", help="Keep hotspot variants.")
    ] = True,
    dry_run: Annotated[
        bool, typer.Option("--dry-run", help="Parse and report counts without writing output.")
    ] = False,
    verbose: VerboseOpt = False,
    quiet: QuietOpt = False,
) -> None:
    """Convert a VEP-annotated SNV/indel VCF into a Curator table."""
    logger = configure_logging(verbose=verbose, quiet=quiet)

    filters = SnvFilterConfig(
        impacts=frozenset(impact) if impact else SnvFilterConfig().impacts,
        keep_filters=frozenset(keep_filter) if keep_filter else None,
        include_splice=include_splice,
        include_hotspot=include_hotspot,
        wgs_min_alt=wgs_min_alt,
    )
    parser = SnvParser(
        vcftype=vcftype,
        oncokb_path=oncokb,
        cgc_path=cgc,
        hotspot_snv_path=hotspot_snv,
        hotspot_indel_path=hotspot_indel,
        curation_path=curation,
        wgs=wgs,
        filters=filters,
    )

    try:
        if dry_run:
            count = sum(1 for _ in parser.parse(input_file))
            logger.info("Dry run: %d rows would be written to %s", count, output)
            typer.echo(f"[dry-run] {count} rows")
            return
        writer = CuratorWriter(parser.columns())
        count = writer.write(parser.parse(input_file), output)
        typer.echo(f"Wrote {count} rows to {output}")
    except Vcf2CuratorError as exc:
        fail(exc)
