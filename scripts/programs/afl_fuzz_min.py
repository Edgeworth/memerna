# Copyright 2022 Eliot Courtney.
from pathlib import Path
from typing import Any

import click
import cloup
from scripts.build.afl_fuzz import AflFuzzCfg
from scripts.build.args import build_cfg_from_args
from scripts.build.args import build_cfg_options
from scripts.util.command import run_shell
from scripts.util.util import fn_args


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
