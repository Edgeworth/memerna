# Copyright 2022 Eliot Courtney.
from pathlib import Path
from typing import Any

import cloup
from rnapy.analysis.fold_accuracy.plotter import FoldAccuracyPlotter
from rnapy.analysis.fold_accuracy.runner import FoldAccuracyRunner
from rnapy.bridge.args import bridge_options
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
def run_fold_accuracy(
    memevault_path: Path,
    dataset: str,
    output_dir: Path,
    memerna: MemeRna,
    **_kwargs: Any,
) -> None:
    memevault = MemeVault(memevault_path, dataset)
    analyser = FoldAccuracyRunner(
        memevault,
        output_dir,
        memerna,
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
def plot_fold_accuracy(input_dir: Path, output_dir: Path) -> None:
    plotter = FoldAccuracyPlotter(input_dir, output_dir)
    plotter.run()
