# Copyright 2022 Eliot Courtney.
from pathlib import Path
from typing import Any

import click
import cloup

from rnapy.build.args import build_cfg_from_args, build_cfg_options
from rnapy.util.command import run_shell
from rnapy.util.util import fn_args


@cloup.command()
@build_cfg_options
@cloup.option("--regenerate/--no-regenerate", default=False)
@cloup.option("--build/--no-build", default=True)
@cloup.option("--test/--no-test", default=False)
@cloup.option("--bench/--no-bench", default=False)
@cloup.option("--bench-output", type=cloup.Path(dir_okay=False, resolve_path=True, path_type=Path))
@cloup.argument("targets", nargs=-1)
def build(
    regenerate: bool,
    build: bool,
    test: bool,
    bench: bool,
    bench_output: Path | None,
    targets: list[str],
    **_kwargs: Any,
) -> None:
    cfg = build_cfg_from_args(**fn_args())

    build_path = cfg.build_path()
    click.echo(build_path)

    cfg.build(targets, build=build, regenerate=regenerate)

    if test:
        run_shell("./run_tests", cwd=build_path)

    if bench:
        cmd = "./run_benchmark"
        if bench_output:
            # Move old benchmark output to .old file, if it exists
            if bench_output.exists():
                bench_output.rename(bench_output.with_suffix(bench_output.suffix + ".old"))
            cmd += f" --benchmark_out={bench_output} --benchmark_out_format=json"
        run_shell(cmd, cwd=build_path)
