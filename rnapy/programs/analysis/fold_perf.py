# Copyright 2022 Eliot Courtney.
from pathlib import Path
from typing import Any

import cloup

from rnapy.analysis.fold_perf import FoldPerfRunner
from rnapy.bridge.args import bridge_options
from rnapy.bridge.linearfold import LinearFold
from rnapy.bridge.memerna01 import MemeRna01
from rnapy.bridge.rnastructure import RNAstructure
from rnapy.bridge.sparsemfefold import SparseMFEFold
from rnapy.bridge.sparsernafold import SparseRNAFolD
from rnapy.bridge.viennarna import ViennaRna
from rnapy.data.args import memevault_options
from rnapy.data.memevault import MemeVault
from rnapy.util.args import init_package_limits, limit_options
from rnapy.util.util import fn_args


@cloup.command()
@bridge_options
@memevault_options
@limit_options
@cloup.option("--dataset", default="random", type=str)
@cloup.option("--num-tries", default=5, type=int)
@cloup.option(
    "--output-path",
    type=cloup.Path(dir_okay=False, file_okay=True, exists=False, path_type=Path),
    required=True,
)
def run_fold_perf(
    num_tries: int,
    memevault_path: Path,
    dataset: str,
    output_path: Path,
    memerna01: MemeRna01,
    linearfold: LinearFold,
    rnastructure: RNAstructure,
    viennarna: ViennaRna,
    sparsemfefold: SparseMFEFold,
    sparsernafold: SparseRNAFolD,
    **_kwargs: Any,
) -> None:
    init_package_limits(**fn_args())
    memevault = MemeVault(memevault_path, dataset)
    analyser = FoldPerfRunner(
        num_tries=num_tries,
        memevault=memevault,
        output_path=output_path,
        memerna01=memerna01,
        linearfold=linearfold,
        rnastructure=rnastructure,
        viennarna=viennarna,
        sparsemfefold=sparsemfefold,
        sparsernafold=sparsernafold,
    )
    analyser.run()
