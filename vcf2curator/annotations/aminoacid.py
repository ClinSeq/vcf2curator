"""Amino-acid notation helpers.

Ported from ``tri_to_single_aa`` in the legacy ``generateIGVnavInput.py``.
"""

from __future__ import annotations

import re

_THREE_TO_ONE: dict[str, str] = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
}

_TOKEN_RE = re.compile(r"([0-9]+|[a-zA-Z\s\n.]+)")


def tri_to_single_aa(hgvsp: str) -> str:
    """Convert three-letter amino-acid codes in an HGVSp string to one-letter.

    Splits the string into alphabetic/numeric tokens, translates pure
    three-letter codes, and preserves any frameshift (``fs``) suffix. Matches
    the legacy behaviour exactly.

    Args:
        hgvsp: The protein-change portion of an HGVSp annotation (without the
            leading ``p.``).

    Returns:
        The HGVSp string with three-letter codes collapsed to one-letter codes.
    """
    res = ""
    tokens = [tok for tok in _TOKEN_RE.split(hgvsp) if tok]

    for token in tokens:
        code = token
        fs_parts: list[str] = []
        # Only translate tokens that are purely alphanumeric (skip '*', '_', etc.).
        if not any(not char.isalnum() for char in code):
            if "fs" in code:
                fs_parts = code.split("fs")
                code = fs_parts[0]
            if not code.isdigit() and code.upper() in _THREE_TO_ONE:
                code = _THREE_TO_ONE[code.upper()]
        res += code
        if fs_parts:
            res += "fs" + fs_parts[1]
    return res
