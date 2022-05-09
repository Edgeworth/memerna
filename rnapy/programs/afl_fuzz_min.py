# Copyright 2022 E.
from pathlib import Path
from typing import Any

import click
import cloup
from rnapy.build.afl_fuzz import AflFuzzCfg
from rnapy.build.args import build_cfg_from_args
from rnapy.build.args import build_cfg_options
from rnapy.util.command import run_shell
from rnapy.util.util import fn_args


@cloup.command()
@build_cfg_options
@cloup.argument(
    "path",
    type=cloup.Path(exists=True, dir_okay=False, resolve_path=True, path_type=Path),
)
def afl_fuzz_min(
    path: Path,
    **_kwargs: Any,
) -> None:
    build_cfg = build_cfg_from_args(**fn_args())
    cfg = AflFuzzCfg(build_cfg=build_cfg)
    cfg.build()

    cmd = cfg.afl_tmin_cmd(path)
    click.echo(f"Running minimisation {cmd}")
    run_shell(cmd, cwd=cfg.bin_path())