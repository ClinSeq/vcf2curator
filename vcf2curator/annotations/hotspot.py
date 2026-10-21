"""Cancer hotspot annotation.

Ported from the hotspot helpers in the legacy ``generateIGVnavInput.py``.
A variant is tagged ``Hotspot`` or ``Warmspot`` based on cancer hotspot
reference tables, plus a small set of hardcoded non-coding hotspot loci.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from vcf2curator.annotations.aminoacid import tri_to_single_aa
from vcf2curator.exceptions import MissingAnnotationError

# Hardcoded non-coding hotspot loci (e.g. TERT promoter) from the legacy script.
_NONCODING_LOCI: frozenset[str] = frozenset({"5:1295228", "5:1295250"})
_NONCODING_INTERVALS: dict[str, tuple[int, int]] = {"7": (116411709, 116411902)}


def load_hotspot(filepath: str | Path, vtype: str) -> pd.DataFrame:
    """Load a cancer hotspot CSV into a DataFrame.

    For SNVs the amino-acid position is expanded into a +/-1 window
    (``start_aa`` / ``end_aa``) used for warmspot matching.

    Args:
        filepath: Path to the hotspot CSV.
        vtype: Either ``"snv"`` or ``"indel"``.

    Returns:
        The hotspot lookup as a DataFrame.

    Raises:
        MissingAnnotationError: If the file cannot be read.
    """
    path = Path(filepath)
    try:
        lookup = pd.read_csv(path)
    except OSError as exc:
        raise MissingAnnotationError(f"Could not read hotspot file {path}: {exc}") from exc

    if vtype == "snv":
        lookup["start_aa"] = lookup["amino_acid_position"] - 1
        lookup["end_aa"] = lookup["amino_acid_position"] + 1
    return lookup


def _hotspot_aapos(gene: str, one_aa: str, aa_range: list[str], lookup: pd.DataFrame) -> str:
    """Classify a single amino-acid position as Hotspot/Warmspot/none."""
    exact = lookup.query(
        "gene == @gene and hgvsp == @one_aa",
        local_dict={"gene": gene, "one_aa": one_aa},
    )
    if len(exact) > 0:
        return "Hotspot"
    window = lookup.query(
        "gene == @gene and start_aa <= @pos <= end_aa",
        local_dict={"gene": gene, "pos": int(aa_range[0])},
    )
    return "Warmspot" if len(window) > 0 else ""


def _hotspot_aarange(gene: str, aa_range: list[str], lookup: pd.DataFrame) -> str:
    """Classify an amino-acid range (inframe indel) as Hotspot/Warmspot/none."""
    scope = {"gene": gene, "start": int(aa_range[0]), "end": int(aa_range[1])}
    both = lookup.query(
        "gene == @gene and (start_aa <= @start <= end_aa) and (start_aa <= @end <= end_aa)",
        local_dict=scope,
    )
    if len(both) > 0:
        return "Hotspot"
    either = lookup.query(
        "gene == @gene and ((start_aa <= @start <= end_aa) or (start_aa <= @end <= end_aa))",
        local_dict=scope,
    )
    return "Warmspot" if len(either) > 0 else ""


def _hotspot_splice(gene: str, aa_range: list[str], gpos: str, lookup: pd.DataFrame) -> str:
    """Classify a splice variant against genomic-position hotspots."""
    res = lookup.query(
        "gene == @gene and amino_acid_position == @pos",
        local_dict={"gene": gene, "pos": int(aa_range[0])},
    )
    if len(res) > 0:
        genomic_pos = res["genomic_position"].iloc[0]
        return "Hotspot" if gpos in genomic_pos else ""
    return ""


def check_noncoding_hotspot(chrom: str, pos: int) -> bool:
    """Return True if a genomic position falls in a known non-coding hotspot.

    Args:
        chrom: Chromosome name (without ``chr`` prefix).
        pos: 1-based genomic position.

    Returns:
        True if the locus is a known non-coding hotspot.
    """
    gpos = f"{chrom}:{pos}"
    if gpos in _NONCODING_LOCI:
        return True
    if chrom in _NONCODING_INTERVALS:
        start, end = _NONCODING_INTERVALS[chrom]
        if start <= pos <= end:
            return True
    return False


def annotate_hotspot(
    gene: str,
    hgvsp: str,
    vartype: str,
    consequence: str,
    lookup: pd.DataFrame,
    gpos: str = "",
) -> str:
    """Annotate a variant with cancer hotspot information.

    Args:
        gene: Gene symbol.
        hgvsp: Full HGVSp string (expected to contain ``p.``).
        vartype: ``"snv"`` or ``"indel"``.
        consequence: VEP consequence string.
        lookup: The hotspot DataFrame for this variant type.
        gpos: ``chrom:pos`` string, used for splice hotspots.

    Returns:
        ``"Hotspot"``, ``"Warmspot"``, or ``""``.
    """
    one_aa = tri_to_single_aa(hgvsp.split("p.")[1])
    aa_range = re.findall(r"\d+", one_aa)
    if not aa_range:
        return ""

    ann = ""
    if vartype == "snv":
        ann = _hotspot_aapos(gene, one_aa, aa_range, lookup)

    if "inframe" in consequence:
        if "_" in one_aa and len(aa_range) > 1:
            ann = _hotspot_aarange(gene, aa_range, lookup)
        else:
            ann = _hotspot_aapos(gene, one_aa, aa_range, lookup)

    if "splice" in consequence and vartype == "snv":
        ann = _hotspot_splice(gene, aa_range, gpos, lookup)

    return ann
