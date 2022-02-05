# Copyright 2022 Eliot Courtney.
from dataclasses import dataclass
from pathlib import Path
import resource
import subprocess
import sys
from typing import Optional
import click

from scripts.util.format import human_size


@dataclass
class CmdLimits:
    time_sec: Optional[int] = None  # limit to time in seconds
    mem_bytes: Optional[int] = None  # limit for rss in bytes


@dataclass
class CmdResult:
    maxrss_bytes: int  # Maximum resident set size in bytes of the process.
    user_sec: float  # User time in seconds.
    sys_sec: float  # System time in seconds.
    real_sec: float  # Real time in seconds.
    ret_code: int  # Return code of the process.
    stdout: str
    stderr: str

    def __str__(self):
        return f"{self.real:.2f}s, {human_size(self.maxrss)} "


def try_cmd(
    *cmd,
    input: Optional[str | bytes] = None,
    return_stdout: bool = True,
    stdout_path: Optional[Path] = None,
    cwd: Optional[Path] = None,
    limits: CmdLimits = CmdLimits(),
):
    if isinstance(input, str):
        input = input.encode("utf-8")

    # Uses GNU time.
    cmd = ["/usr/bin/time", "-f", "%e %U %S %M"] + list(cmd)

    def pre_exec():
        if limits.mem_bytes:
            resource.setrlimit(resource.R.AS, (limits.mem_bytes, limits.mem_bytes))
        if limits.time_sec:
            resource.setrlimit(resource.RLIMIT_CPU, (limits.time_sec, limits.time_sec))

    if stdout_path:
        stdout = open(stdout_path, "w")
    else:
        stdout = subprocess.PIPE if return_stdout else subprocess.DEVNULL
    stdin = subprocess.PIPE if input else None
    with subprocess.Popen(
        cmd,
        shell=False,
        stdin=stdin,
        stdout=stdout,
        stderr=subprocess.PIPE,
        cwd=cwd,
        preexec_fn=pre_exec,
    ) as proc:
        proc_stdout, proc_stderr = proc.communicate(input=input)
        ret_code = proc.wait()
        if proc_stdout:
            proc_stdout = proc_stdout.decode("UTF-8")
        if proc_stderr:
            proc_stderr = proc_stderr.decode("UTF-8")

        last_line = proc_stderr.strip().split("\n")[-1].split(" ")
        real_sec, user_sec, sys_sec, maxrss_kb = (float(i) for i in last_line)

        if stdout_path:
            stdout.flush()
            stdout.close()

            # We may want to not return the stdout if it's too big.
            if return_stdout:
                proc_stdout = stdout_path.read_text()

        return CmdResult(
            stdout=proc_stdout,
            stderr=proc_stderr,
            ret_code=ret_code,
            real_sec=real_sec,
            user_sec=user_sec,
            sys_sec=sys_sec,
            maxrss_bytes=maxrss_kb * 1024,
        )


def run_cmd(
    *cmd,
    input: Optional[str | bytes] = None,
    return_stdout: bool = True,
    stdout_path: Optional[Path] = None,
    cwd: Optional[Path] = None,
    limits: CmdLimits = CmdLimits(),
):
    res = try_cmd(
        *cmd,
        input=input,
        return_stdout=return_stdout,
        stdout_path=stdout_path,
        cwd=cwd,
        limits=limits,
    )
    if res.ret_code != 0:
        click.echo(f"Running `{cmd}' failed with ret code {res.ret_code}.")
        click.echo(f"stderr: {res.stderr}")
        sys.exit(1)
    return res
