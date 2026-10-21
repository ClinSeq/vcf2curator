"""OncoKB "all annotated variants" lookup.

Ported from ``loadOncoKB`` in the legacy ``generateIGVnavInput.py``.
"""

from __future__ import annotations

from pathlib import Path

from vcf2curator.exceptions import MissingAnnotationError

# OncoKB record: oncogenicity, mutation effect, PMIDs (as originally emitted).
OncoKBEntry = list[str]
OncoKBLookup = dict[str, dict[str, OncoKBEntry]]


def load_oncokb(filepath: str | Path) -> OncoKBLookup:
    """Load the OncoKB all-variants TSV into a gene -> protein_change lookup.

    The file is a tab-delimited export where (0-based) columns are:
    3=gene, 5=protein change, 6=oncogenicity, 7=mutation effect, 8=PMIDs.

    Args:
        filepath: Path to the OncoKB all-annotated-variants TSV.

    Returns:
        Mapping of gene symbol -> {protein_change -> [oncogenicity, effect, pmids]}.

    Raises:
        MissingAnnotationError: If the file cannot be read.
    """
    path = Path(filepath)
    lookup: OncoKBLookup = {}
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                data = line.rstrip("\n").split("\t")
                if len(data) < 8:
                    continue
                gene = data[3]
                protein_change = data[5]
                oncogenicity = data[6]
                mutation_effect = data[7]
                pmids = data[8] if len(data) > 8 else ""
                lookup.setdefault(gene, {})[protein_change] = [
                    oncogenicity,
                    mutation_effect,
                    pmids,
                ]
    except OSError as exc:
        raise MissingAnnotationError(f"Could not read OncoKB file {path}: {exc}") from exc
    return lookup


def lookup_oncogenicity(
    oncokb: OncoKBLookup,
    gene: str,
    amino_acids: list[str],
    protein_positions: list[str],
) -> str | OncoKBEntry:
    """Return the OncoKB entry for a variant, or an empty string if none.

    Reproduces the legacy protein-change key of ``ref_aa + position + alt_aa``.

    Args:
        oncokb: The loaded OncoKB lookup.
        gene: Gene symbol.
        amino_acids: VEP ``Amino_acids`` split on ``/``.
        protein_positions: VEP ``Protein_position`` split on ``/``.

    Returns:
        The ``[oncogenicity, effect, pmids]`` entry, or ``""`` when not found.
    """
    if gene not in oncokb:
        return ""
    if len(amino_acids) > 1 and len(protein_positions) > 0:
        protein_change = amino_acids[0] + protein_positions[0] + amino_acids[1]
        entry = oncokb[gene].get(protein_change)
        if entry is not None:
            return entry
    return ""
