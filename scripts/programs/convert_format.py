# Copyright 2016 E.
from typing import Any

import click
import cloup
from scripts.data.args import data_options
from scripts.data.args import rna_from_args
from scripts.model.parse.rna_parser import RnaParser
from scripts.util.util import fn_args


@cloup.command(aliases=["conv"])
@data_options
@cloup.option("-k", "--kind", type=cloup.Choice(["db", "ct"]), default="db")
def convert_format(kind: str, **_kwargs: Any) -> None:
    rna = rna_from_args(**fn_args())
    if kind == "ct":
        click.echo(RnaParser.to_ct_file(rna))
    if kind == "db":
        click.echo(RnaParser.to_db_file(rna))
