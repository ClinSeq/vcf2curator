"""GRIDSS breakend parsing and simple SV type classification.

Replaces the legacy ``gridss_svannotate.R`` (``simpleEventType`` from
StructuralVariantAnnotation). Each breakend record's ALT already describes its
partner, so an event is classified from a single record without looking up
the mate.

Strand follows StructuralVariantAnnotation: ``+`` means the breakend keeps the
sequence before the position, ``-`` the sequence after it.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import StrEnum

# ALT breakend notations (VCF 4.2 section 5.4): t[p[, t]p], ]p]t, [p[t.
_ALT_PATTERNS: tuple[tuple[re.Pattern[str], str, str, bool], ...] = (
    # (pattern, local strand, partner strand, local bases come first)
    (re.compile(r"^([A-Za-z]+)\[([^:\[\]]+):(\d+)\[$"), "+", "-", True),
    (re.compile(r"^([A-Za-z]+)\]([^:\[\]]+):(\d+)\]$"), "+", "+", True),
    (re.compile(r"^\]([^:\[\]]+):(\d+)\]([A-Za-z]+)$"), "-", "+", False),
    (re.compile(r"^\[([^:\[\]]+):(\d+)\[([A-Za-z]+)$"), "-", "-", False),
)

# Fraction of |svLen| the inserted sequence must reach to call an insertion.
_INS_FRACTION = 0.7


class SvType(StrEnum):
    """Simple structural variant event types."""

    DEL = "DEL"
    DUP = "DUP"
    INS = "INS"
    INV = "INV"
    TRA = "TRA"


@dataclass(frozen=True)
class Breakend:
    """One side of a breakpoint, with its partner parsed from ALT.

    Attributes:
        chrom: Chromosome of this breakend.
        pos: 1-based position of this breakend.
        strand: Orientation of this breakend (``+`` or ``-``).
        partner_chrom: Chromosome of the partner breakend.
        partner_pos: 1-based position of the partner breakend.
        partner_strand: Orientation of the partner breakend.
        inserted: Untemplated sequence inserted at the breakpoint.
    """

    chrom: str
    pos: int
    strand: str
    partner_chrom: str
    partner_pos: int
    partner_strand: str
    inserted: str


def parse_breakend(chrom: str, pos: int, ref: str, alt: str) -> Breakend | None:
    """Parse a paired breakend ALT.

    Args:
        chrom: Record chromosome.
        pos: Record 1-based position.
        ref: Record REF allele.
        alt: Record ALT allele.

    Returns:
        The parsed breakend, or ``None`` for single breakends (``A.``/``.A``)
        and any ALT that is not paired breakend notation.
    """
    for pattern, strand, partner_strand, local_first in _ALT_PATTERNS:
        match = pattern.match(alt)
        if match is None:
            continue
        if local_first:
            bases, partner_chrom, partner_pos = match.groups()
            inserted = bases[len(ref) :]
        else:
            partner_chrom, partner_pos, bases = match.groups()
            inserted = bases[: len(bases) - len(ref)]
        return Breakend(
            chrom=chrom,
            pos=pos,
            strand=strand,
            partner_chrom=partner_chrom,
            partner_pos=int(partner_pos),
            partner_strand=partner_strand,
            inserted=inserted,
        )
    return None


def _is_deletion_like(bnd: Breakend) -> bool:
    """Return True when the breakends join across removed reference sequence."""
    return (bnd.pos < bnd.partner_pos) != (bnd.strand == "-")


def sv_length(bnd: Breakend) -> int | None:
    """Return the StructuralVariantAnnotation-style svLen, or None across chromosomes.

    Deletion-like events are negative; inserted bases are added on top.
    """
    if bnd.chrom != bnd.partner_chrom:
        return None
    distance = abs(bnd.pos - bnd.partner_pos) - 1
    signed = -distance if _is_deletion_like(bnd) else distance
    return signed + len(bnd.inserted)


def simple_type(bnd: Breakend) -> SvType:
    """Classify a breakend into a simple event type.

    Same decision order as ``simpleEventType`` in the legacy R script.
    """
    if bnd.chrom != bnd.partner_chrom:
        return SvType.TRA
    if bnd.strand == bnd.partner_strand:
        return SvType.INV
    length = sv_length(bnd)
    if length is not None and len(bnd.inserted) >= abs(length) * _INS_FRACTION:
        return SvType.INS
    if _is_deletion_like(bnd):
        return SvType.DEL
    return SvType.DUP
