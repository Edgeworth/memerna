# Copyright 2022 Eliot Courtney.
import os
import shutil
from typing import Any

import click
import cloup
from scripts.build.args import build_cfg_from_args
from scripts.build.args import build_cfg_options
from scripts.build.cmake import generate_cmake
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

    if regenerate and build_path.exists():
        shutil.rmtree(build_path)

    click.echo(build_path)
    if not os.path.exists(build_path):
        regenerate = True

    if regenerate:
        click.echo("Regenerating cmake files.")
        generate_cmake(cfg)

    if build:
        # Add stack protector etc to catch non-crashing memory bugs.
        run_shell(f"AFL_HARDEN=1 make -j$(($(nproc)-1)) {' '.join(targets)}", cwd=cfg.build_path())

    if test:
        run_shell("./run_tests", cwd=cfg.build_path())
