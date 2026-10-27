"""Typer CLI entry point: builds the app and registers each command."""

from __future__ import annotations

from typing import Annotated

import typer

from vcf2curator import __version__
from vcf2curator.commands.igv_relocate import igv_relocate
from vcf2curator.commands.igv_session import igv_session
from vcf2curator.commands.snvs import snvs
from vcf2curator.commands.svs import svs
from vcf2curator.commands.validate import validate

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Transform raw bioinformatics pipeline output into Curator-ready tables.",
)

app.command()(validate)
app.command()(snvs)
app.command()(svs)
app.command(name="igv-session")(igv_session)
app.command(name="igv-relocate")(igv_relocate)


def version_callback(value: bool) -> None:
    """Print the version and exit."""
    if value:
        typer.echo(f"vcf2curator {__version__}")
        raise typer.Exit()


@app.callback()
def _main(
    version: Annotated[
        bool,
        typer.Option(
            "--version", callback=version_callback, is_eager=True, help="Show version and exit."
        ),
    ] = False,
) -> None:
    """vcf2curator command-line interface."""


def main() -> None:
    """Console-script entry point."""
    app()


if __name__ == "__main__":
    main()
