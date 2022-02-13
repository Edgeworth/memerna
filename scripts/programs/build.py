# Copyright 2022 Eliot Courtney.
from typing import Any

import click
import cloup
from scripts.build.args import build_cfg_from_args
from scripts.build.args import build_cfg_options
from scripts.util.command import run_shell
from scripts.util.util import fn_args


@cloup.command()
@build_cfg_options
@cloup.option("--regenerate/--no-regenerate", default=False)
@cloup.option("--build/--no-build", default=True)
@cloup.option("--test/--no-test", default=False)
@cloup.argument("targets", nargs=-1)
def build(
    regenerate: bool,
    build: bool,
    test: bool,
    targets: list[str],
    **_kwargs: Any,
) -> None:
    cfg = build_cfg_from_args(**fn_args())

    build_path = cfg.build_path()
    click.echo(build_path)

    if build:
        cfg.build(targets, regenerate=regenerate)

    if test:
        run_shell("./run_tests", cwd=build_path)
