# Copyright 2022 E.
from pathlib import Path
from typing import Any

import cloup

from rnapy.analysis.subopt_perf.plotter import SuboptPerfPlotter
from rnapy.analysis.subopt_perf.runner import SuboptPerfRunner
from rnapy.bridge.args import bridge_options
from rnapy.bridge.memerna import MemeRna
from rnapy.bridge.rnastructure import RNAstructure
from rnapy.bridge.viennarna import ViennaRna
from rnapy.data.args import memevault_options
from rnapy.data.memevault import MemeVault


@cloup.command()
@bridge_options
@memevault_options
@cloup.option("--dataset", default="random", type=str)
@cloup.option("--time-sec-limit", type=int, required=False)
@cloup.option("--mem-bytes-limit", type=int, required=False)
@cloup.option("--num-tries", default=5, type=int)
@cloup.option(
    "--output-dir",
    type=cloup.Path(dir_okay=True, file_okay=False, exists=True, path_type=Path),
    required=True,
)
@cloup.option("--delta/--no-delta", default=True)
def run_subopt_perf(
    time_sec_limit: int | None,
    mem_bytes_limit: int | None,
    num_tries: int,
    memevault_path: Path,
    dataset: str,
    output_dir: Path,
    delta: bool,
    memerna: MemeRna,
    rnastructure: RNAstructure,
    viennarna: ViennaRna,
    **_kwargs: Any,
) -> None:
    memevault = MemeVault(memevault_path, dataset)
    analyser = SuboptPerfRunner(
        time_sec_limit=time_sec_limit,
        mem_bytes_limit=mem_bytes_limit,
        num_tries=num_tries,
        memevault=memevault,
        output_dir=output_dir,
        delta=delta,
        memerna=memerna,
        rnastructure=rnastructure,
        viennarna=viennarna,
    )
    analyser.run()


@cloup.command()
@cloup.option(
    "--input-dir",
    type=cloup.Path(dir_okay=True, file_okay=False, exists=True, path_type=Path),
    required=True,
)
@cloup.option(
    "--output-dir",
    type=cloup.Path(dir_okay=True, file_okay=False, exists=True, path_type=Path),
    required=True,
)
def plot_subopt_perf(input_dir: Path, output_dir: Path) -> None:
    plotter = SuboptPerfPlotter(input_dir, output_dir)
    plotter.run()
