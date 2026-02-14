# Copyright 2022 Eliot Courtney.
import os
import signal
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
    compress_stdout: bool = False,
    cwd: Path | None = None,
    extra_env: dict[str, str] | None = None,
    limits: CmdLimits | None = None,
) -> CmdResult:
    limits = limits or CmdLimits()
    if isinstance(stdin_inp, str):
        stdin_inp = stdin_inp.encode("utf-8")

    cmd_list: list[str] = list(cmd)
    if limits.cpu_affinity is not None:
        cmd_list = [
            "taskset",
            "-c",
            ",".join(str(c) for c in sorted(limits.cpu_affinity)),
            *cmd_list,
        ]
    prlimit_args: list[str] = []
    if limits.time_sec is not None:
        prlimit_args.append(f"--cpu={limits.time_sec}")
    if limits.mem_bytes is not None:
        prlimit_args.append(f"--as={limits.mem_bytes}")
    if prlimit_args:
        cmd_list = ["prlimit", *prlimit_args, "--", *cmd_list]

    # Uses GNU time.
    cmd = ("/usr/bin/time", "-f", "%e %U %S %M", *cmd_list)

    env = os.environ.copy()
    if extra_env is not None:
        env.update(extra_env)

    stdout: Any
    zstd_proc: subprocess.Popen[bytes] | None = None
    stdout_file: Any = None
    if stdout_path is not None:
        if compress_stdout:
            stdout_file = stdout_path.open("wb")
            zstd_proc = subprocess.Popen(
                ["zstd", "--fast", "-q"],  # noqa: S607
                stdin=subprocess.PIPE,
                stdout=stdout_file,
            )
            stdout = zstd_proc.stdin
        else:
            stdout = stdout_path.open("wb")
    else:
        stdout = subprocess.PIPE if stdout_to_str else subprocess.DEVNULL
    stdin = subprocess.PIPE if stdin_inp is not None else None

    CMD_STR_LIM = 500
    cmd_str = " ".join(cmd)
    if len(cmd_str) > CMD_STR_LIM:
        cmd_str = cmd_str[: CMD_STR_LIM // 2] + "..." + cmd_str[-CMD_STR_LIM // 2 :]
    print(
        f"try_cmd: {cmd_str}, cwd: {cwd}, extra_env: {extra_env}, "
        f"stdout_to_str: {stdout_to_str}, stdout_path: {stdout_path}, "
        f"compress_stdout: {compress_stdout}, limits: {limits}"
    )
    with subprocess.Popen(
        cmd, shell=False, stdin=stdin, stdout=stdout, stderr=subprocess.PIPE, cwd=cwd, env=env
    ) as proc:
        # Close parent's copy of zstd's stdin so zstd sees EOF when main proc exits.
        if zstd_proc is not None and zstd_proc.stdin is not None:
            zstd_proc.stdin.close()

        stdout_bytes, stderr_bytes = proc.communicate(input=stdin_inp)
        ret_code = proc.wait()
        stdout_str = stdout_bytes.decode("utf-8") if stdout_bytes else ""
        stderr_str = stderr_bytes.decode("utf-8") if stderr_bytes else ""

        last_line = stderr_str.strip().rsplit("\n", maxsplit=1)[-1].split(" ")
        real_sec, user_sec, sys_sec, maxrss_kb = (float(i) for i in last_line)

        if stdout_path is not None:
            if zstd_proc is not None:
                zstd_proc.wait()
                stdout_file.close()
                if zstd_proc.returncode != 0:
                    raise RuntimeError(
                        f"zstd compression failed with return code {zstd_proc.returncode}"
                    )
            else:
                stdout.flush()
                stdout.close()

            # We may want to not return the stdout if it's too big.
            if stdout_to_str:
                if compress_stdout:
                    decomp = subprocess.run(
                        ["zstd", "-dc", str(stdout_path)],  # noqa: S607
                        capture_output=True,
                        text=True,
                        check=False,
                    )
                    if decomp.returncode != 0:
                        raise RuntimeError(f"zstd -dc failed: {decomp.stderr}")
                    stdout_str = decomp.stdout
                else:
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
    compress_stdout: bool = False,
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
        compress_stdout=compress_stdout,
        cwd=cwd,
        extra_env=extra_env,
        limits=limits,
    )
    if res.ret_code != 0:
        # If killed by SIGINT/SIGTERM (e.g. Ctrl-C), propagate as KeyboardInterrupt
        # so callers don't treat it as a failure.
        if res.ret_code < 0 and -res.ret_code in (signal.SIGINT, signal.SIGTERM):
            raise KeyboardInterrupt
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
