"""``validate`` command: check a VCF is usable before converting it."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from vcf2curator.commands._shared import QuietOpt, VerboseOpt, fail
from vcf2curator.commands.snvs import SnvParser
from vcf2curator.curator.schema import VcfType
from vcf2curator.exceptions import Vcf2CuratorError
from vcf2curator.logging_config import configure_logging


def validate(
    input_file: Annotated[
        Path,
        typer.Option(
            "--input", "-i", exists=True, dir_okay=False, help="VEP-annotated VCF to validate."
        ),
    ],
    verbose: VerboseOpt = False,
    quiet: QuietOpt = False,
) -> None:
    """Check a VCF is readable and VEP-annotated with the required CSQ fields."""
    logger = configure_logging(verbose=verbose, quiet=quiet)
    try:
        # Validation only inspects the VCF header, so no annotation resources are needed.
        SnvParser(vcftype=VcfType.SOMATIC, oncokb_path="").validate(input_file)
    except Vcf2CuratorError as exc:
        fail(exc)
    logger.info("%s is a valid VEP-annotated VCF", input_file)
    typer.echo("OK")
