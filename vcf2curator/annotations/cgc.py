"""Cancer Gene Census (CGC) lookup.

Ported from ``loadCGC`` in the legacy ``generateIGVnavInput.py``. The SNV script
keys the lookup by Ensembl gene ID; the SV command keys it by gene symbol.
"""

from __future__ import annotations

from pathlib import Path

from vcf2curator.exceptions import MissingAnnotationError

# CGC entry: (annotation, gene_symbol), keyed by Ensembl gene ID.
CgcLookup = dict[str, tuple[str, str]]


def load_cgc(filepath: str | Path) -> CgcLookup:
    """Load the Cancer Gene Census TSV keyed by Ensembl gene ID.

    Column layout (0-based, header skipped): 0=gene symbol, 1=Ensembl ID,
    4=annotation.

    Args:
        filepath: Path to the CGC TSV.

    Returns:
        Mapping of Ensembl gene ID -> (annotation, gene symbol).

    Raises:
        MissingAnnotationError: If the file cannot be read.
    """
    path = Path(filepath)
    lookup: CgcLookup = {}
    try:
        with path.open("r", encoding="utf-8") as handle:
            handle.readline()  # header
            for line in handle:
                data = line.rstrip("\n").split("\t")
                if len(data) < 5:
                    continue
                gene_symbol = data[0]
                ensembl_id = data[1]
                annotation = data[4]
                lookup[ensembl_id] = (annotation, gene_symbol)
    except OSError as exc:
        raise MissingAnnotationError(f"Could not read CGC file {path}: {exc}") from exc
    return lookup


def load_cgc_by_symbol(filepath: str | Path) -> dict[str, str]:
    """Load the Cancer Gene Census TSV keyed by gene symbol (used by ``svs``).

    Args:
        filepath: Path to the CGC TSV.

    Returns:
        Mapping of gene symbol -> annotation.

    Raises:
        MissingAnnotationError: If the file cannot be read.
    """
    path = Path(filepath)
    lookup: dict[str, str] = {}
    try:
        with path.open("r", encoding="utf-8") as handle:
            handle.readline()  # header
            for line in handle:
                data = line.rstrip("\n").split("\t")
                if len(data) < 5:
                    continue
                lookup[data[0]] = data[4]
    except OSError as exc:
        raise MissingAnnotationError(f"Could not read CGC file {path}: {exc}") from exc
    return lookup
