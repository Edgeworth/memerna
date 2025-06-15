# Copyright 2022 Eliot Courtney.
import multiprocessing
import shlex
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


def run_fuzz(cfg: AflFuzzCfg, window: libtmux.Window) -> None:
    cfg.build()
    cmd = cfg.afl_fuzz_cmd()
    cwd = cfg.bin_path()
    click.echo(f"Running fuzz {cmd} in {cwd}")

    # Run in given tmux window:
    pane = window.attached_pane
    if not pane:
        raise RuntimeError(f"Window {window} has no attached pane")
    pane.send_keys(f"cd {shlex.quote(str(cwd))}", enter=True, suppress_history=True)
    pane.send_keys(cmd, enter=True, suppress_history=True)


@cloup.command()
@build_cfg_options
@afl_fuzz_cfg_options
@cloup.option(
    "--num-procs",
    default=multiprocessing.cpu_count() - 2,
    help="Number of fuzzing configurations to run.",
)
def afl_fuzz(num_procs: int, **_kwargs: Any) -> None:
    build_cfg = build_cfg_from_args(**fn_args())
    afl_cfg = build_afl_fuzz_cfg_from_args(build_cfg, **fn_args())
    cfgs = afl_fuzz_cfgs(afl_cfg, num_procs)

    # Use libtmux to set up a session with a random name.
    session = None
    try:
        server = libtmux.Server()
        session_name = "afl_fuzz"
        session = server.new_session(session_name, kill_session=True)
        windows = []
        for i in range(len(cfgs)):
            window = session.new_window(attach=False, window_name=f"window_{i}")
            windows.append(window)

        with multiprocessing.Pool(len(cfgs)) as pool:
            pool.starmap(run_fuzz, zip(cfgs, windows, strict=True))

        click.echo("Attaching session")
        session.attach_session()
    except Exception:
        click.echo("Error occurred, killing session")
        if session is not None:
            session.kill_session()
        raise
