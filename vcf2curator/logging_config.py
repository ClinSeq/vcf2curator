"""Consistent logging setup shared by every subcommand."""

from __future__ import annotations

import logging

LOGGER_NAME = "vcf2curator"

_LOG_FORMAT = "%(levelname)s %(asctime)s %(name)s - %(message)s"


def configure_logging(*, verbose: bool = False, quiet: bool = False) -> logging.Logger:
    """Configure and return the package logger.

    Args:
        verbose: If True, set the level to DEBUG.
        quiet: If True, set the level to ERROR. Takes precedence over ``verbose``.

    Returns:
        The configured package logger.
    """
    if quiet:
        level = logging.ERROR
    elif verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logging.basicConfig(level=level, format=_LOG_FORMAT)
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(level)
    return logger


def get_logger() -> logging.Logger:
    """Return the package logger without reconfiguring handlers."""
    return logging.getLogger(LOGGER_NAME)
