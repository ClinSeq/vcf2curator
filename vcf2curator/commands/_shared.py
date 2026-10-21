"""Option types and error handling shared by every CLI command."""

from __future__ import annotations

from typing import Annotated, NoReturn

import typer

from vcf2curator.exceptions import (
    InputValidationError,
    MissingAnnotationError,
    SchemaError,
    Vcf2CuratorError,
)
from vcf2curator.logging_config import configure_logging

# Exit codes distinguish failure modes for scripting.
_EXIT_CODES: dict[type[Vcf2CuratorError], int] = {
    InputValidationError: 3,
    MissingAnnotationError: 4,
    SchemaError: 5,
}

VerboseOpt = Annotated[bool, typer.Option("--verbose", "-v", help="Enable debug logging.")]
QuietOpt = Annotated[bool, typer.Option("--quiet", "-q", help="Only log errors.")]


def fail(exc: Vcf2CuratorError) -> NoReturn:
    """Log a domain error and exit with its dedicated exit code."""
    logger = configure_logging()
    logger.error("%s", exc)
    raise typer.Exit(code=_EXIT_CODES.get(type(exc), 1))
