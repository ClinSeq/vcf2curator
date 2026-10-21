"""Supporting-read collection, the SAME_START_READS filter, and evidence BAMs.

Replaces the legacy ``generate_evidence_bam.py``. Instead of building an
in-memory name index over the whole BAM (``pysam.IndexedReads``), the BAM is
streamed once and only alignments of the wanted read names are kept.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Iterable
from pathlib import Path

import pysam

from vcf2curator.exceptions import InputValidationError
from vcf2curator.logging_config import get_logger

logger = get_logger()

# Read name -> its non-supplementary alignments, in BAM order.
AlignmentsByName = dict[str, list[pysam.AlignedSegment]]


def collect_alignments(bam_path: str | Path, read_names: set[str]) -> AlignmentsByName:
    """Return the non-supplementary alignments of ``read_names`` from one BAM pass.

    Args:
        bam_path: Path to the BAM the variants were called from.
        read_names: Read names to collect.

    Returns:
        Mapping of read name -> alignments. Names absent from the BAM are omitted.

    Raises:
        InputValidationError: If the BAM cannot be opened.
    """
    try:
        bam = pysam.AlignmentFile(str(bam_path), "rb")
    except (OSError, ValueError) as exc:
        raise InputValidationError(f"Could not open BAM {bam_path}: {exc}") from exc

    found: AlignmentsByName = defaultdict(list)
    with bam:
        for read in bam.fetch(until_eof=True):
            if not read.is_supplementary and read.query_name in read_names:
                found[read.query_name].append(read)
    logger.info("Found %d of %d supporting reads in %s", len(found), len(read_names), bam_path)
    return dict(found)


def _signature(alignments: list[pysam.AlignedSegment]) -> str:
    """Summarize where a read's alignments sit, e.g. ``1:100:250-1:400:550``."""
    return "-".join(f"{a.reference_name}:{a.reference_start}:{a.reference_end}" for a in alignments)


def has_same_start_reads(
    read_names: Iterable[str],
    alignments: AlignmentsByName,
    *,
    fraction: float,
    min_distinct: int,
) -> bool:
    """Return True if an event's supporting reads look like a single artefact.

    An event fails when one alignment signature accounts for at least
    ``fraction`` of its read names, or when fewer than ``min_distinct``
    distinct signatures are seen.

    Args:
        read_names: The event's supporting read names (INFO/BPNAMES).
        alignments: Alignments collected by :func:`collect_alignments`.
        fraction: Share of reads one signature may reach before failing.
        min_distinct: Minimum number of distinct signatures required.
    """
    names = set(read_names)
    signatures = Counter(_signature(alignments[n]) for n in names if n in alignments)
    cutoff = len(names) * fraction
    return len(signatures) < min_distinct or any(c >= cutoff for c in signatures.values())


def write_evidence_bam(
    template_bam: str | Path,
    alignments: AlignmentsByName,
    read_names: Iterable[str],
    output: str | Path,
) -> int:
    """Write the alignments of ``read_names`` to a sorted, indexed BAM.

    Args:
        template_bam: BAM whose header is copied to the output.
        alignments: Alignments collected by :func:`collect_alignments`.
        read_names: Read names to include.
        output: Destination BAM path; ``<output>.bai`` is written next to it.

    Returns:
        The number of alignments written.
    """
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    unsorted = output.with_name(output.name + ".unsorted.bam")
    count = 0
    with pysam.AlignmentFile(str(template_bam), "rb") as template:
        with pysam.AlignmentFile(str(unsorted), "wb", template=template) as out:
            for name in sorted(set(read_names)):
                for read in alignments.get(name, []):
                    out.write(read)
                    count += 1
    pysam.sort("-o", str(output), str(unsorted))
    unsorted.unlink()
    pysam.index(str(output))
    logger.info("Wrote %d alignments to evidence BAM %s", count, output)
    return count
