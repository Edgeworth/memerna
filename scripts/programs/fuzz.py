# Copyright 2022 Eliot Courtney.
import multiprocessing

from typing import Any

import click
import cloup
from scripts.build.args import build_cfg_from_args
from scripts.build.args import build_cfg_options
from scripts.build.cmake import generate_cmake

from scripts.build.fuzz import FuzzCfg, build_fuzz_cfgs
from scripts.util.command import run_shell
from scripts.util.util import fn_args


def run_fuzz(cfg: FuzzCfg) -> None:
    cmd = cfg.afl_fuzz_cmd()
    click.echo(f"Running fuzz {cmd}")
    run_shell(cmd, cwd=cfg.bin_path())


@cloup.command()
@build_cfg_options
@cloup.option(
    "--num-procs",
    default=multiprocessing.cpu_count() - 1,
    help="Number of fuzzing configurations to run.",
)
def fuzz(
    num_procs: int,
    **_kwargs: Any,
) -> None:
    build_cfg = build_cfg_from_args(**fn_args())
    cfgs = build_fuzz_cfgs(build_cfg, num_procs)

    generate_cmake(build_cfg)
    for cfg in cfgs:
        cfg.build()

    with multiprocessing.Pool(len(cfgs)) as pool:
        pool.map(run_fuzz, cfgs)
