# Copyright 2022 Eliot Courtney.
from dataclasses import dataclass
import os
from pathlib import Path
import resource
import subprocess
from typing import Any

import click
from scripts.util.format import human_size


@dataclass
class CmdLimits:
    time_sec: int | None = None  # limit to time in seconds
    mem_bytes: int | None = None  # limit for rss in bytes


@dataclass
class CmdResult:
    maxrss_bytes: int  # Maximum resident set size in bytes of the process.
    user_sec: float  # User time in seconds.
    sys_sec: float  # System time in seconds.
    real_sec: float  # Real time in seconds.
    ret_code: int  # Return code of the process.
    stdout: str
    stderr: str

    def __str__(self) -> str:
        return f"{self.real_sec:.2f}s, {human_size(self.maxrss_bytes)} "


def try_cmd(  # pylint: disable=too-many-locals
    *cmd: str,
    inp: str | bytes | None = None,
    return_stdout: bool = True,
    stdout_path: Path | None = None,
    cwd: Path | None = None,
    extra_env: dict[str, str] | None = None,
    limits: CmdLimits = CmdLimits(),
) -> CmdResult:
    if isinstance(inp, str):
        inp = inp.encode("utf-8")

    # Uses GNU time.
    cmd = ("/usr/bin/time", "-f", "%e %U %S %M") + cmd

    def preexec_fn() -> None:
        if limits.mem_bytes:
            resource.setrlimit(resource.RLIMIT_AS, (limits.mem_bytes, limits.mem_bytes))
        if limits.time_sec:
            resource.setrlimit(resource.RLIMIT_CPU, (limits.time_sec, limits.time_sec))

    env = os.environ.copy()
    if extra_env:
        env.update(extra_env)

    stdout: Any
    if stdout_path:
        stdout = open(stdout_path, "wb")  # pylint: disable=consider-using-with
    else:
        stdout = subprocess.PIPE if return_stdout else subprocess.DEVNULL
    stdin = subprocess.PIPE if inp else None
    with subprocess.Popen(  # pylint: disable=subprocess-popen-preexec-fn
        cmd,
        shell=False,
        stdin=stdin,
        stdout=stdout,
        stderr=subprocess.PIPE,
        cwd=cwd,
        env=env,
        preexec_fn=preexec_fn,
    ) as proc:
        stdout_bytes, stderr_bytes = proc.communicate(input=inp)
        ret_code = proc.wait()
        stdout_str = stdout_bytes.decode("utf-8") if stdout_bytes else ""
        stderr_str = stderr_bytes.decode("utf-8") if stderr_bytes else ""

        last_line = stderr_str.strip().rsplit("\n", maxsplit=1)[-1].split(" ")
        real_sec, user_sec, sys_sec, maxrss_kb = (float(i) for i in last_line)

        if stdout_path:
            stdout.flush()
            stdout.close()

            # We may want to not return the stdout if it's too big.
            if return_stdout:
                stdout_str = stdout_path.read_text()

        return CmdResult(
            stdout=stdout_str,
            stderr=stderr_str,
            ret_code=ret_code,
            real_sec=real_sec,
            user_sec=user_sec,
            sys_sec=sys_sec,
            maxrss_bytes=int(maxrss_kb) * 1024,
        )


def run_cmd(
    *cmd: str,
    inp: str | bytes | None = None,
    return_stdout: bool = True,
    stdout_path: Path | None = None,
    cwd: Path | None = None,
    extra_env: dict[str, str] | None = None,
    limits: CmdLimits = CmdLimits(),
) -> CmdResult:
    res = try_cmd(
        *cmd,
        inp=inp,
        return_stdout=return_stdout,
        stdout_path=stdout_path,
        cwd=cwd,
        extra_env=extra_env,
        limits=limits,
    )
    if res.ret_code != 0:
        click.echo(f"Running `{cmd}' failed with ret code {res.ret_code}.")
        click.echo(f"stderr: {res.stderr}")
        raise RuntimeError(f"Shell command failed: {cmd}")
    return res


def run_shell(cmd: str, cwd: Path | None = None) -> None:
    prev_cwd = os.getcwd()
    if cwd is not None:
        os.chdir(cwd)
    ret_code = os.system(cmd)
    os.chdir(prev_cwd)
    if ret_code != 0:
        click.echo(f"Running `{cmd}' failed with ret code {ret_code}.")
        raise RuntimeError(f"Shell command failed: {cmd}")
