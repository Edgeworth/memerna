# Copyright 2016 Eliot Courtney.
import click
import cloup
from scripts.data.args import data_options
from scripts.data.args import rna_from_args
from scripts.util.util import fn_args


@cloup.command(aliases=["conv"])
@data_options
@cloup.option("-t", "--type", type=cloup.Choice(["db", "ct"]), default="db")
def convert_format(type, **kwargs):
    rna = rna_from_args(**fn_args())
    match type:
        case "ct":
            click.echo(rna.to_ct_file())
        case "db":
            click.echo(rna.to_db_file())
