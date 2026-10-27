"""``igv-relocate`` command: rewrite path prefixes in existing IGV session files.

Rewrites ``Resource/@path`` and ``Track/@id`` values in place. A prefix matches
a whole path or a leading directory: ``/data/a`` matches ``/data/a/x.bam`` but
not ``/data/ab/x.bam``. The first matching ``--replace`` wins.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Annotated

import typer

from vcf2curator.commands._shared import QuietOpt, VerboseOpt, fail
from vcf2curator.exceptions import InputValidationError, Vcf2CuratorError
from vcf2curator.logging_config import configure_logging, get_logger

logger = get_logger()

# Element tag -> attribute holding a file path.
_PATH_ATTRIBUTES = {"Resource": "path", "Track": "id"}

PrefixMap = list[tuple[str, str]]


def parse_replacements(values: list[str]) -> PrefixMap:
    """Parse ``OLD=NEW`` values into ``(old, new)`` prefix pairs.

    Raises:
        InputValidationError: If a value has no ``=`` or an empty OLD prefix.
    """
    pairs: PrefixMap = []
    for value in values:
        old, sep, new = value.partition("=")
        old = old.rstrip("/")
        if not sep or not old:
            raise InputValidationError(f"--replace expects OLD=NEW with a non-empty OLD: {value!r}")
        pairs.append((old, new.rstrip("/")))
    return pairs


def relocate_value(value: str, replacements: PrefixMap) -> str:
    """Return ``value`` with the first matching prefix replaced."""
    for old, new in replacements:
        if value == old or value.startswith(old + "/"):
            return new + value[len(old) :]
    return value


def relocate_session(root: ET.Element, replacements: PrefixMap) -> int:
    """Rewrite path attributes under ``root`` and return the number changed."""
    changed = 0
    for element in root.iter():
        attribute = _PATH_ATTRIBUTES.get(element.tag)
        value = element.get(attribute) if attribute else None
        if attribute is None or value is None:
            continue
        new_value = relocate_value(value, replacements)
        if new_value != value:
            element.set(attribute, new_value)
            changed += 1
    return changed


def _parse_session(path: Path) -> ET.ElementTree[ET.Element]:
    try:
        return ET.parse(path)
    except (OSError, ET.ParseError) as exc:
        raise InputValidationError(f"Could not read IGV session {path}: {exc}") from exc


def igv_relocate(
    sessions: Annotated[
        list[Path],
        typer.Argument(exists=True, dir_okay=False, help="IGV session XML files to rewrite."),
    ],
    replace: Annotated[
        list[str],
        typer.Option("--replace", help="Path prefix mapping OLD=NEW (repeatable)."),
    ],
    dry_run: Annotated[
        bool, typer.Option("--dry-run", help="Report changes without writing files.")
    ] = False,
    verbose: VerboseOpt = False,
    quiet: QuietOpt = False,
) -> None:
    """Rewrite file path prefixes in IGV session files, e.g. after moving results."""
    configure_logging(verbose=verbose, quiet=quiet)
    try:
        replacements = parse_replacements(replace)
        trees = {path: _parse_session(path) for path in sessions}
    except Vcf2CuratorError as exc:
        fail(exc)

    for path, tree in trees.items():
        changed = relocate_session(tree.getroot(), replacements)
        if changed == 0:
            logger.info("%s: no paths matched", path)
            continue
        if dry_run:
            logger.info("Dry run: would rewrite %d paths in %s", changed, path)
            continue
        tree.write(path, encoding="UTF-8", xml_declaration=True)
        logger.info("Rewrote %d paths in %s", changed, path)
