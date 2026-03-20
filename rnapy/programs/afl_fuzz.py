# Copyright 2022 Eliot Courtney.
import multiprocessing
import shlex
from collections import defaultdict
from typing import Any

import click
import cloup
import libtmux

from rnapy.build.afl_fuzz import AflFuzzCfg, afl_fuzz_cfgs
from rnapy.build.args import (
    afl_fuzz_cfg_options,
    build_afl_fuzz_cfg_from_args,
    build_cfg_from_args,
    build_cfg_options,
)
from rnapy.util.util import fn_args


def build_fuzz(cfg: AflFuzzCfg) -> None:
    cfg.build()


def launch_fuzz(cfg: AflFuzzCfg, window: libtmux.Window) -> None:
    cmd = cfg.afl_fuzz_cmd()
    cwd = cfg.bin_path()
    click.echo(f"Running fuzz {cmd} in {cwd}")

    # Run in given tmux window:
    pane = window.active_pane
    if not pane:
        raise RuntimeError(f"Window {window} has no attached pane")
    pane.send_keys(f"cd {shlex.quote(str(cwd))}", enter=True, suppress_history=True)
    pane.send_keys(cmd, enter=True, suppress_history=True)


@cloup.command()
@build_cfg_options
@afl_fuzz_cfg_options
@cloup.option(
    "--num-procs",
    default=max(1, multiprocessing.cpu_count() - 2),
    help="Number of fuzzing configurations to run.",
)
def afl_fuzz(num_procs: int, **_kwargs: Any) -> None:
    build_cfg = build_cfg_from_args(**fn_args())
    afl_cfg = build_afl_fuzz_cfg_from_args(build_cfg, **fn_args())
    cfgs = afl_fuzz_cfgs(afl_cfg, num_procs)

    # Group configs by cmake build directory.
    by_build_path: dict[str, list[AflFuzzCfg]] = defaultdict(list)
    for cfg in cfgs:
        by_build_path[str(cfg.build_cfg.build_path())].append(cfg)

    # Build unique configs in parallel.
    representatives = [group[0] for group in by_build_path.values()]
    with multiprocessing.Pool(len(representatives)) as pool:
        pool.map(build_fuzz, representatives)

    # Copy binaries for remaining same-kind fuzzers that didn't build.
    for group in by_build_path.values():
        for cfg in group[1:]:
            cfg.build()

    # Launch each fuzzer in its own tmux window.
    session = None
    try:
        server = libtmux.Server()
        session_name = "afl_fuzz"
        session = server.new_session(session_name, kill_session=True)
        windows = []
        for i in range(len(cfgs)):
            window = session.new_window(attach=False, window_name=f"window_{i}")
            windows.append(window)

        for cfg, window in zip(cfgs, windows, strict=True):
            launch_fuzz(cfg, window)

        click.echo("Attaching session")
        session.attach()
    except Exception:
        click.echo("Error occurred, killing session")
        if session is not None:
            session.kill()
        raise
