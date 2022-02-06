# Copyright 2022 Eliot Courtney.
from inspect import signature
import inspect
from pathlib import Path
import cloup
from scripts.bridge.args import bridge_options
from scripts.data.args import data_options, rna_from_args
from scripts.model.args import (
    energy_cfg_from_args,
    energy_options,
    subopt_options,
    subopt_cfg_from_args,
)
from scripts.model.config import EnergyCfg, SuboptCfg
from scripts.model.rna import Rna
from scripts.util.util import fn_args


def run_fold(program, rna):
    frna, res = program.fold(rna)
    print(f"Folding {rna.name} with {program}: {frna.db()}\n  {res}")


def run_efn(program, rna):
    energy, res = program.efn(rna)
    print(f"Energy of {rna.name} with {program}: {energy:f}\n  {res}")


def run_suboptimal(program, rna, delta, subopt_max_print):
    subopts, res = program.suboptimal(rna, delta, None)
    if res.ret:
        print("Execution failed")
        return
    print(f"{len(subopts)} suboptimal structures of {rna.name} with {program} - {res}")
    subopt_subset = subopts
    if subopt_max_print > 0:
        subopt_subset = subopts[:subopt_max_print]
    for energy, structure in subopt_subset:
        print(f"{energy:f} {structure.db()}")


def run_partition(program, rna):
    pass  # TODO: implement


@cloup.command()
@energy_options
@subopt_options
@bridge_options
@data_options
@cloup.option("-e", "--efn", type=bool, default=False, help="Run efn the given RNA")
@cloup.option("-f", "--fold", type=bool, default=False, help="Fold the given RNA")
@cloup.option("-p", "--partition", type=bool, default=False, help="Run partition on the given RNA")
@cloup.option(
    "-s", "--subopt", type=bool, default=False, help="Run suboptimal folding on the given RNA"
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
        ]
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
    programs = []
    for program in programs:
        programs.append(kwargs[program])

    for program in programs:
        if fold:
            run_fold(program, rna)
        if efn:
            run_efn(program, rna)
        if partition:
            run_partition(program, rna)
        if subopt:
            run_suboptimal(program, rna, energy_cfg, subopt_cfg)
