"""Cancer-specific curation tag lookup.

Ported from ``annotate_curate_tag`` in the legacy ``generateIGVnavInput.py``.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from vcf2curator.exceptions import MissingAnnotationError


def load_curation(filepath: str | Path) -> pd.DataFrame:
    """Load the curation annotation CSV.

    Args:
        filepath: Path to the curation CSV (columns include ``prefix``, ``gene``,
            ``type`` and a comment/tag column).

    Returns:
        The curation table as a DataFrame.

    Raises:
        MissingAnnotationError: If the file cannot be read.
    """
    path = Path(filepath)
    try:
        return pd.read_csv(path)
    except OSError as exc:
        raise MissingAnnotationError(f"Could not read curation file {path}: {exc}") from exc


def annotate_curate_tag(
    curation: pd.DataFrame,
    prefix: str,
    gene: str,
    vartype: str,
) -> str:
    """Return the curation tag for a variant, or an empty string.

    Args:
        curation: The loaded curation DataFrame.
        prefix: Project prefix (derived from the input filename).
        gene: Gene symbol.
        vartype: ``"somatic"`` or ``"germline"``.

    Returns:
        ``"<gene> <tag>"`` when a matching curation row exists, else ``""``.
    """
    res = curation.query("prefix == @prefix and gene == @gene and type == @vartype")
    if len(res) > 0:
        curate_tag = res.iloc[0, 3]
        return " ".join([gene, str(curate_tag)])
    return ""
