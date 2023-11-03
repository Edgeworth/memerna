# Copyright 2022 Eliot Courtney.
from typing import Any

import click
import cloup

from rnapy.bridge.args import bridge_options
from rnapy.bridge.rnapackage import RnaPackage
from rnapy.data.args import data_options, rna_from_args
from rnapy.model.args import (
    energy_cfg_from_args,
    energy_options,
    subopt_cfg_from_args,
    subopt_options,
)
from rnapy.util.util import fn_args


@cloup.command()
@energy_options
@subopt_options
@bridge_options
@data_options
@cloup.option("-e", "--efn/--no-efn", default=False, help="Run efn the given RNA")
@cloup.option("-f", "--fold/--no-fold", default=False, help="Fold the given RNA")
@cloup.option(
    "-p", "--partition/--no-partition", default=False, help="Run partition on the given RNA"
)
@cloup.option(
    "-s", "--subopt/--no-subopt", default=False, help="Run suboptimal folding on the given RNA"
)
@cloup.option(
    "-p",
    "--programs",
    multiple=True,
    type=cloup.Choice(["memerna", "rnastructure", "sparsemfefold", "viennarna"]),
    help="Programs to run",
)
def harness(  # pylint: disable=too-many-locals
    programs: list[str], efn: bool, fold: bool, partition: bool, subopt: bool, **kwargs: Any
) -> None:
    energy_cfg = energy_cfg_from_args(**fn_args())
    subopt_cfg = subopt_cfg_from_args(**fn_args())
    rna = rna_from_args(**fn_args())
    packages: list[RnaPackage] = []

    for program in programs:
        package = kwargs.get(program)
        if package is None:
            raise click.UsageError(f"Path for program {program} not found")
        packages.append(kwargs[program])

    for p in packages:
        if efn:
            energy, res = p.efn(rna, energy_cfg)
            click.echo(f"Energy of RNA {rna.name} with {p}: {energy:f}\n  {res}")
        if fold:
            frna, res = p.fold(rna, energy_cfg)
            click.echo(f"Folding RNA {rna.name} with {p}: {frna.db()}\n  {res}")
        if partition:
            raise NotImplementedError
        if subopt:
            subopts, res = p.subopt(rna, energy_cfg, subopt_cfg)
            click.echo(f"{len(subopts)} suboptimal structures of RNA {rna.name} with {p} - {res}")
            for rna in subopts:
                click.echo(f"{rna.energy} {rna.db()}")
