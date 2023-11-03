# Copyright 2022 Eliot Courtney.
from pathlib import Path
from typing import Any

import click
import cloup

from rnapy.build.args import (
    afl_fuzz_cfg_options,
    build_afl_fuzz_cfg_from_args,
    build_cfg_from_args,
    build_cfg_options,
)
from rnapy.util.command import run_shell
from rnapy.util.util import fn_args


@cloup.command()
@build_cfg_options
@afl_fuzz_cfg_options
@cloup.argument(
    "paths",
    type=cloup.Path(exists=True, dir_okay=False, resolve_path=True, path_type=Path),
    nargs=-1,
)
def afl_fuzz_min(paths: list[Path], **_kwargs: Any) -> None:
    build_cfg = build_cfg_from_args(**fn_args())
    afl_cfg = build_afl_fuzz_cfg_from_args(build_cfg, **fn_args())
    afl_cfg.build()

    for path in paths:
        cmd = afl_cfg.afl_tmin_cmd(path)
        click.echo(f"Running minimisation {cmd}")
        run_shell(cmd, cwd=afl_cfg.bin_path())
