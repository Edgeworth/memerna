# Copyright 2022 Eliot Courtney.
import click
import cloup
from scripts.bridge.args import bridge_options
from scripts.bridge.rnapackage import RnaPackage
from scripts.data.args import data_options
from scripts.data.args import rna_from_args
from scripts.model.args import energy_cfg_from_args
from scripts.model.args import energy_options
from scripts.model.args import subopt_cfg_from_args
from scripts.model.args import subopt_options
from scripts.util.util import fn_args


@cloup.command()
@energy_options
@subopt_options
@bridge_options
@data_options
@cloup.option("-e", "--efn", type=bool, default=False, help="Run efn the given RNA")
@cloup.option("-f", "--fold", type=bool, default=False, help="Fold the given RNA")
@cloup.option("-p", "--partition", type=bool, default=False, help="Run partition on the given RNA")
@cloup.option(
    "-s",
    "--subopt",
    type=bool,
    default=False,
    help="Run suboptimal folding on the given RNA",
)
@cloup.option(
    "-p",
    "--programs",
    multiple=True,
    type=cloup.Choice(
        [
            "memerna",
            "rnastructure",
            "sparsemfefold",
            "viennarna",
        ],
    ),
    help="Programs to run",
)
def harness(
    programs: list[str],
    efn: bool,
    fold: bool,
    partition: bool,
    subopt: bool,
    **kwargs,
):
    energy_cfg = energy_cfg_from_args(**fn_args())
    subopt_cfg = subopt_cfg_from_args(**fn_args())
    rna = rna_from_args(**fn_args())
    programs: list[RnaPackage] = []
    for program in programs:
        programs.append(kwargs[program])

    for program in programs:
        if fold:
            frna, res = program.fold(rna, energy_cfg)
            click.echo(f"Folding {rna.name} with {program}: {frna.db()}\n  {res}")
        if efn:
            energy, res = program.efn(rna, energy_cfg)
            click.echo(f"Energy of {rna.name} with {program}: {energy:f}\n  {res}")
        if partition:
            raise NotImplementedError
        if subopt:
            subopts, res = program.subopt(rna, energy_cfg, subopt_cfg)
            click.echo(f"{len(subopts)} suboptimal structures of {rna.name} with {program} - {res}")
            for energy, structure in subopts:
                click.echo(f"{energy} {structure.db()}")
