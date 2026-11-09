"""Tests for rewriting path prefixes in IGV session files."""

from __future__ import annotations

import xml.etree.ElementTree as ET

import pytest

from vcf2curator.commands.igv_relocate import parse_replacements, relocate_session, relocate_value
from vcf2curator.exceptions import InputValidationError

_SESSION = """<Session genome="hg19">
    <Resources>
        <Resource path="/old/data/t.bam"/>
        <Resource path="/old/database/x.bed"/>
        <Resource name="Ensembl Genes" path="https://example.org/ensGene.txt.gz"/>
    </Resources>
    <Panel name="t.bam">
        <Track id="/old/data/t.bam_coverage" name="t.bam Coverage"/>
        <Track id="Reference sequence" name="Reference sequence"/>
    </Panel>
</Session>"""


def test_parse_replacements_strips_trailing_slash() -> None:
    assert parse_replacements(["/old/=/new/", "/a=/b"]) == [("/old", "/new"), ("/a", "/b")]


@pytest.mark.parametrize("value", ["/old", "=/new", "/"])
def test_parse_replacements_rejects_bad_values(value: str) -> None:
    with pytest.raises(InputValidationError, match="OLD=NEW"):
        parse_replacements([value])


def test_relocate_value_matches_directory_boundary() -> None:
    mapping = [("/old/data", "/new")]
    assert relocate_value("/old/data/t.bam", mapping) == "/new/t.bam"
    assert relocate_value("/old/data", mapping) == "/new"
    assert relocate_value("/old/database/x.bed", mapping) == "/old/database/x.bed"


def test_relocate_value_first_match_wins() -> None:
    mapping = [("/old/data", "/first"), ("/old", "/second")]
    assert relocate_value("/old/data/t.bam", mapping) == "/first/t.bam"


def test_relocate_session_rewrites_resources_and_track_ids() -> None:
    root = ET.fromstring(_SESSION)
    assert relocate_session(root, [("/old/data", "/mnt/data")]) == 2
    assert [r.attrib["path"] for r in root.iter("Resource")] == [
        "/mnt/data/t.bam",
        "/old/database/x.bed",
        "https://example.org/ensGene.txt.gz",
    ]
    assert [t.attrib["id"] for t in root.iter("Track")] == [
        "/mnt/data/t.bam_coverage",
        "Reference sequence",
    ]
