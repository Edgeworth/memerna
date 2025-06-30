# Copyright 2022 Eliot Courtney.
from pathlib import Path
from typing import Any

import cloup

from rnapy.analysis.fold_accuracy import FoldAccuracyRunner
from rnapy.bridge.args import bridge_options
from rnapy.bridge.memerna import MemeRna
from rnapy.data.args import memevault_options
from rnapy.data.memevault import MemeVault
from rnapy.util.args import init_package_limits, limit_options
from rnapy.util.util import fn_args


@cloup.command()
@bridge_options
@memevault_options
@limit_options
@cloup.option("--dataset", default="random", type=str)
@cloup.option(
    "--output-path", type=cloup.Path(file_okay=True, dir_okay=False, path_type=Path), required=True
)
def run_fold_accuracy(
    memevault_path: Path, dataset: str, output_path: Path, memerna: MemeRna, **_kwargs: Any
) -> None:
    init_package_limits(**fn_args())
    memevault = MemeVault(memevault_path, dataset)
    analyser = FoldAccuracyRunner(memevault=memevault, output_path=output_path, memerna=memerna)
    analyser.run()
