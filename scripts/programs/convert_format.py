# Copyright 2016 Eliot Courtney.
from typing import Any
import click
import cloup
from scripts.data.args import data_options
from scripts.data.args import rna_from_args
from scripts.util.util import fn_args


@cloup.command(aliases=["conv"])
@data_options
@cloup.option("-k", "--kind", type=cloup.Choice(["db", "ct"]), default="db")
def convert_format(kind: str, **_kwargs: Any) -> None:
    rna = rna_from_args(**fn_args())
    if kind == "ct":
        click.echo(rna.to_ct_file())
    if kind == "db":
        click.echo(rna.to_db_file())
