# Copyright 2022 Eliot Courtney.
import subprocess
from typing import Any

import cloup

from rnapy.build.afl_fuzz import afl_fuzz_cfg_by_index
from rnapy.build.args import (
    afl_fuzz_cfg_options,
    afl_fuzz_index_option,
    build_afl_fuzz_cfg_from_args,
    build_cfg_from_args,
    build_cfg_options,
)
from rnapy.util.util import fn_args


@cloup.command()
@build_cfg_options
@afl_fuzz_cfg_options
@afl_fuzz_index_option
@cloup.argument("testcase", type=str)
def afl_fuzz_run(testcase: str, index: int, **_kwargs: Any) -> None:
    build_cfg = build_cfg_from_args(**fn_args())
    afl_cfg = build_afl_fuzz_cfg_from_args(build_cfg, **fn_args())
    afl_cfg = afl_fuzz_cfg_by_index(afl_cfg, index)
    afl_cfg.build()
    res = subprocess.run(
        afl_cfg.fuzz_argv(), input=testcase + "\n", cwd=afl_cfg.bin_path(), text=True, check=False
    )
    if res.returncode != 0:
        raise RuntimeError(f"Fuzz target exited with code {res.returncode}.")
