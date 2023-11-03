# Copyright 2016 Eliot Courtney.
from typing import Any

import click
import cloup

from rnapy.data.args import data_options, rna_from_args
from rnapy.model.parse.rna_parser import RnaParser
from rnapy.util.util import fn_args


@cloup.command(aliases=["conv"])
@data_options
@cloup.option("-k", "--kind", type=cloup.Choice(["db", "ct"]), default="db")
def convert_format(kind: str, **_kwargs: Any) -> None:
    rna = rna_from_args(**fn_args())
    if kind == "ct":
        click.echo(RnaParser.to_ct_file(rna))
    if kind == "db":
        click.echo(RnaParser.to_db_file(rna))
