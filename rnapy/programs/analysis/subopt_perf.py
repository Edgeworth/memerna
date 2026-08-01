# Copyright 2022 Eliot Courtney.
from pathlib import Path
from typing import Any

import cloup
from cloup.constraints import require_all

from rnapy.analysis.subopt_perf import SuboptPerfRunner
from rnapy.bridge.args import memerna_path_option, rnastructure_path_option, viennarna_path_option
from rnapy.bridge.memerna import MemeRna
from rnapy.bridge.rnastructure import RNAstructure
from rnapy.bridge.viennarna import ViennaRna
from rnapy.data.args import memevault_path_option
from rnapy.data.memevault import MemeVault
from rnapy.util.args import init_package_limits, limit_options
from rnapy.util.util import fn_args


@cloup.command()
@cloup.option_group(
    "Bridge options",
    memerna_path_option,
    rnastructure_path_option,
    viennarna_path_option,
    constraint=require_all,
)
@memevault_path_option(required=True)
@limit_options
@cloup.option("--dataset", default="random", type=str)
@cloup.option("--num-tries", default=5, type=cloup.IntRange(min=1))
@cloup.option("--rna-lengths", default=None, type=cloup.IntRange(min=1), multiple=True)
@cloup.option(
    "-j", "--jobs", default=1, type=cloup.IntRange(min=1), help="number of parallel workers"
)
@cloup.option(
    "--output-path",
    type=cloup.Path(dir_okay=False, file_okay=True, exists=False, path_type=Path),
    required=True,
)
def run_subopt_perf(
    num_tries: int,
    rna_lengths: tuple[int, ...],
    jobs: int,
    memevault_path: Path,
    dataset: str,
    output_path: Path,
    memerna: MemeRna,
    rnastructure: RNAstructure,
    viennarna: ViennaRna,
    **_kwargs: Any,
) -> None:
    init_package_limits(**fn_args())
    memevault = MemeVault(memevault_path, dataset)
    analyser = SuboptPerfRunner(
        num_tries=num_tries,
        rna_lengths=rna_lengths,
        jobs=jobs,
        memevault=memevault,
        output_path=output_path,
        memerna=memerna,
        rnastructure=rnastructure,
        viennarna=viennarna,
    )
    analyser.run()
