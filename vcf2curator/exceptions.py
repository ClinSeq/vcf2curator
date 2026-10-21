"""Domain-specific exceptions for vcf2curator.

Every user-facing failure should be one of these so the CLI can turn it into a
clear message and a distinct exit code instead of a raw traceback.
"""

from __future__ import annotations


class Vcf2CuratorError(Exception):
    """Base class for all vcf2curator errors."""


class InputValidationError(Vcf2CuratorError):
    """Raised when an input file is missing, unreadable, or malformed."""


class MissingAnnotationError(Vcf2CuratorError):
    """Raised when a required annotation resource is missing or malformed."""


class SchemaError(Vcf2CuratorError):
    """Raised when an output row does not match the Curator schema."""
