"""Curator output schema definitions.

These column lists ARE the Curator input contract: they are taken verbatim from
the headers emitted by the legacy ``generateIGVnavInput.py`` script (somatic and
germline variants) and ``generateIGVnavInput_SV.py`` (structural variants, minus
the TOOL column). Column order and names are significant and must not change
without a corresponding change on the Curator side.
"""

from __future__ import annotations

from enum import StrEnum

# Somatic SNV/indel table (legacy somatic header).
SOMATIC_COLUMNS: tuple[str, ...] = (
    "CHROM",
    "START",
    "END",
    "REF",
    "ALT",
    "CALL",
    "TAG",
    "NOTES",
    "GENE",
    "ENSEMBLID",
    "IMPACT",
    "CONSEQUENCE",
    "TRANSCRIPT",
    "HGVSc",
    "HGVSp",
    "AMINO_ACID",
    "HOTSPOT",
    "T_DP",
    "T_ALT",
    "T_VAF",
    "N_DP",
    "N_ALT",
    "N_VAF",
    "CLIN_SIG",
    "RSID",
    "gnomAD",
    "BRCAEx",
    "OncoKB",
    "CGC_ANN",
    "CURATE",
    "NUM_TOOLS",
)

# Germline SNV/indel table (legacy germline header).
GERMLINE_COLUMNS: tuple[str, ...] = (
    "CHROM",
    "START",
    "END",
    "REF",
    "ALT",
    "CALL",
    "TAG",
    "NOTES",
    "GENE",
    "ENSEMBLID",
    "IMPACT",
    "CONSEQUENCE",
    "TRANSCRIPT",
    "HGVSc",
    "HGVSp",
    "AMINO_ACID",
    "HOTSPOT",
    "N_DP",
    "N_ALT",
    "N_VAF",
    "T_VAF",
    "CLIN_SIG",
    "RSID",
    "gnomAD",
    "BRCAEx",
    "OncoKB",
    "CGC_ANN",
    "CURATE",
)

# Structural variant table (legacy SV summary header without the TOOL column).
SV_COLUMNS: tuple[str, ...] = (
    "CHROM_A",
    "START_A",
    "END_A",
    "CHROM_B",
    "START_B",
    "END_B",
    "IGV_COORD",
    "SVTYPE",
    "SV_LENGTH",
    "SUPPORT_READS",
    "VAF",
    "SDID",
    "SAMPLE",
    "GENE_A",
    "GENE_B",
    "GENE_A-GENE_B-sorted",
    "CGC_ANN",
    "CURATOR",
    "CURATE",
)


class VcfType(StrEnum):
    """The kind of variant call set being curated."""

    SOMATIC = "somatic"
    GERMLINE = "germline"


def columns_for(vcftype: VcfType) -> tuple[str, ...]:
    """Return the Curator column order for the given VCF type.

    Args:
        vcftype: Somatic or germline.

    Returns:
        The ordered tuple of column names for that variant type.
    """
    return SOMATIC_COLUMNS if vcftype is VcfType.SOMATIC else GERMLINE_COLUMNS
