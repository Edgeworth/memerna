# Copyright 2022 E.
import os
import resource
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import click

from rnapy.util.format import human_size


@dataclass
class CmdLimits:
    time_sec: int | None = None  # limit to time in seconds
    mem_bytes: int | None = None  # limit for rss in bytes
    cpu_affinity: set[int] | None = None  # set of CPU IDs to run using


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


def try_cmd(
    *cmd: str,
    stdin_inp: str | bytes | None = None,
    stdout_to_str: bool = True,
    stdout_path: Path | None = None,
    cwd: Path | None = None,
    extra_env: dict[str, str] | None = None,
    limits: CmdLimits | None = None,
) -> CmdResult:
    limits = limits or CmdLimits()
    if isinstance(stdin_inp, str):
        stdin_inp = stdin_inp.encode("utf-8")

    # Uses GNU time.
    cmd = ("/usr/bin/time", "-f", "%e %U %S %M", *cmd)

    def preexec_fn() -> None:
        if limits.mem_bytes is not None:
            resource.setrlimit(resource.RLIMIT_AS, (limits.mem_bytes, limits.mem_bytes))
        if limits.time_sec is not None:
            resource.setrlimit(resource.RLIMIT_CPU, (limits.time_sec, limits.time_sec))
        if limits.cpu_affinity is not None:
            os.sched_setaffinity(0, limits.cpu_affinity)

    env = os.environ.copy()
    if extra_env is not None:
        env.update(extra_env)

    stdout: Any
    if stdout_path is not None:
        stdout = stdout_path.open("wb")
    else:
        stdout = subprocess.PIPE if stdout_to_str else subprocess.DEVNULL
    stdin = subprocess.PIPE if stdin_inp else None

    CMD_STR_LIM = 500
    cmd_str = " ".join(cmd)
    if len(cmd_str) > CMD_STR_LIM:
        cmd_str = cmd_str[: CMD_STR_LIM // 2] + "..." + cmd_str[-CMD_STR_LIM // 2 :]
    print(
        f"try_cmd: {cmd_str}, cwd: {cwd}, extra_env: {extra_env}, "
        f"stdout_to_str: {stdout_to_str}, stdout_path: {stdout_path}, limits: {limits}"
    )
    with subprocess.Popen(
        cmd,
        shell=False,
        stdin=stdin,
        stdout=stdout,
        stderr=subprocess.PIPE,
        cwd=cwd,
        env=env,
        preexec_fn=preexec_fn,  # noqa: PLW1509
    ) as proc:
        stdout_bytes, stderr_bytes = proc.communicate(input=stdin_inp)
        ret_code = proc.wait()
        stdout_str = stdout_bytes.decode("utf-8") if stdout_bytes else ""
        stderr_str = stderr_bytes.decode("utf-8") if stderr_bytes else ""

        last_line = stderr_str.strip().rsplit("\n", maxsplit=1)[-1].split(" ")
        real_sec, user_sec, sys_sec, maxrss_kb = (float(i) for i in last_line)

        if stdout_path is not None:
            stdout.flush()
            stdout.close()

            # We may want to not return the stdout if it's too big.
            if stdout_to_str:
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
    stdin_inp: str | bytes | None = None,
    stdout_to_str: bool = True,
    stdout_path: Path | None = None,
    cwd: Path | None = None,
    extra_env: dict[str, str] | None = None,
    limits: CmdLimits | None = None,
) -> CmdResult:
    limits = limits or CmdLimits()
    res = try_cmd(
        *cmd,
        stdin_inp=stdin_inp,
        stdout_to_str=stdout_to_str,
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


def run_shell(cmd: str, cwd: Path | None = None, extra_env: dict[str, str] | None = None) -> None:
    prev_cwd = Path.cwd()
    if cwd is not None:
        os.chdir(cwd)
    if extra_env is not None:
        cmd = f"{' '.join(f'{k}={v}' for k, v in extra_env.items())} {cmd}"
    print(f"cmd: {cmd}, cwd: {cwd}")
    ret_code = os.system(cmd)
    os.chdir(prev_cwd)
    if ret_code != 0:
        click.echo(f"Running `{cmd}' failed with ret code {ret_code}.")
        raise RuntimeError(f"Shell command failed: {cmd}")
