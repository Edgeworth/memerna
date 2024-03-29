# Copyright 2022 E.
from pathlib import Path
from typing import Any

import cloup

from rnapy.analysis.fold_perf.plotter import FoldPerfPlotter
from rnapy.analysis.fold_perf.runner import FoldPerfRunner
from rnapy.bridge.args import bridge_options
from rnapy.bridge.linearfold import LinearFold
from rnapy.bridge.memerna import MemeRna
from rnapy.bridge.rnastructure import RNAstructure
from rnapy.bridge.sparsemfefold import SparseMfeFold
from rnapy.bridge.viennarna import ViennaRna
from rnapy.data.args import memevault_options
from rnapy.data.memevault import MemeVault


@cloup.command()
@bridge_options
@memevault_options
@cloup.option("--dataset", default="random", type=str)
@cloup.option(
    "--output-dir",
    type=cloup.Path(dir_okay=True, file_okay=False, exists=True, path_type=Path),
    required=True,
)
def run_fold_perf(
    memevault_path: Path,
    dataset: str,
    output_dir: Path,
    memerna: MemeRna,
    linearfold: LinearFold,
    rnastructure: RNAstructure,
    viennarna: ViennaRna,
    sparsemfefold: SparseMfeFold,
    **_kwargs: Any,
) -> None:
    memevault = MemeVault(memevault_path, dataset)
    analyser = FoldPerfRunner(
        memevault, output_dir, memerna, linearfold, rnastructure, viennarna, sparsemfefold
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
def plot_fold_perf(input_dir: Path, output_dir: Path) -> None:
    plotter = FoldPerfPlotter(input_dir, output_dir)
    plotter.run()
