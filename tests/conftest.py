"""Shared pytest fixtures."""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def fixtures_dir() -> Path:
    """Return the directory holding sample input files."""
    return FIXTURES


# Read-pair start positions per SV event; "dup" reads all start together so the
# event fails SAME_START_READS.
_EVENT_READ_STARTS = {
    "del": (900, 950, 1010, 1050),
    "dup": (2950, 2950, 2950, 2950),
    "inv": (4900, 4950, 5010, 5050),
    "ins": (6900, 6950, 7010, 7050),
    "tra": (7900, 7950, 8010, 8050),
}


def _read(header: pysam.AlignmentHeader, name: str, start: int, flag: int) -> pysam.AlignedSegment:
    read = pysam.AlignedSegment(header)
    read.query_name = name
    read.flag = flag
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = 60
    read.cigarstring = "50M"
    read.query_sequence = "A" * 50
    read.query_qualities = pysam.qualitystring_to_array("I" * 50)
    return read


@pytest.fixture
def sv_bam(tmp_path: Path) -> Path:
    """Build a small sorted, indexed BAM holding the GRIDSS fixture's supporting reads."""
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": 20000}]}
    )
    reads = []
    for event, starts in _EVENT_READ_STARTS.items():
        for i, start in enumerate(starts, start=1):
            name = f"{event}_r{i}"
            reads.append(_read(header, name, start, flag=0x1 | 0x40))
            reads.append(_read(header, name, start + 300, flag=0x1 | 0x80))
    # A supplementary alignment and an unrelated read, both of which must be ignored.
    reads.append(_read(header, "del_r1", 15000, flag=0x1 | 0x40 | 0x800))
    reads.append(_read(header, "unrelated", 500, flag=0))

    unsorted = tmp_path / "unsorted.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for read in reads:
            out.write(read)
    bam = tmp_path / "tumor.bam"
    pysam.sort("-o", str(bam), str(unsorted))
    pysam.index(str(bam))
    return bam
