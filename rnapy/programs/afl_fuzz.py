# Copyright 2022 Eliot Courtney.
import multiprocessing
from typing import Any

import click
import cloup
from rnapy.build.afl_fuzz import afl_fuzz_cfgs
from rnapy.build.afl_fuzz import AflFuzzCfg
from rnapy.build.args import build_cfg_from_args
from rnapy.build.args import build_cfg_options
from rnapy.util.command import run_shell
from rnapy.util.util import fn_args


def run_fuzz(cfg: AflFuzzCfg) -> None:
    cmd = cfg.afl_fuzz_cmd()
    click.echo(f"Running fuzz {cmd}")
    run_shell(cmd, cwd=cfg.bin_path())


@cloup.command()
@build_cfg_options
@cloup.option(
    "--num-procs",
    default=multiprocessing.cpu_count() - 2,
    help="Number of fuzzing configurations to run.",
)
@cloup.option(
    "--max-len",
    default=-1,
    help="Max size of sequences to fuzz",
)
def afl_fuzz(
    num_procs: int,
    max_len: int,
    **_kwargs: Any,
) -> None:
    build_cfg = build_cfg_from_args(**fn_args())
    cfgs = afl_fuzz_cfgs(build_cfg, num_procs, max_len)

    for cfg in cfgs:
        cfg.build()

    with multiprocessing.Pool(len(cfgs)) as pool:
        pool.map(run_fuzz, cfgs)
