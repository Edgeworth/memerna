# Copyright 2022 Eliot Courtney.
import functools
import os
import shutil
import signal
import subprocess
import sys
import tempfile
from contextlib import ExitStack, suppress
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Any, override

import click

from rnapy.util.format import human_size


@dataclass
class CmdLimits:
    time_sec: int | None = None  # limit to time in seconds
    mem_bytes: int | None = None  # limit for process memory in bytes
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

    @override
    def __str__(self) -> str:
        return f"{self.real_sec:.2f}s, {human_size(self.maxrss_bytes)} "


@functools.cache
def _can_use_systemd_user_scope() -> bool:
    if sys.platform != "linux":
        return False
    if shutil.which("systemd-run") is None:
        return False
    try:
        res = subprocess.run(
            ["systemd-run", "--user", "--scope", "--quiet", "true"],  # noqa: S607
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
            timeout=5,
        )
    except OSError, subprocess.SubprocessError:
        return False
    return res.returncode == 0


@functools.cache
def _warn_systemd_user_scope_unavailable() -> None:
    click.echo(
        "Warning: systemd-run --user --scope is unavailable; falling back to prlimit --as, "
        "which limits virtual address space rather than resident memory.",
        err=True,
    )


def _cmd_list_with_limits(cmd: tuple[str, ...], limits: CmdLimits) -> list[str]:
    cmd_list: list[str] = list(cmd)
    if limits.cpu_affinity is not None:
        cmd_list = [
            "taskset",
            "-c",
            ",".join(str(c) for c in sorted(limits.cpu_affinity)),
            *cmd_list,
        ]

    use_systemd_memory_limit = limits.mem_bytes is not None and _can_use_systemd_user_scope()

    prlimit_args: list[str] = []
    if limits.time_sec is not None:
        prlimit_args.append(f"--cpu={limits.time_sec}")
    if limits.mem_bytes is not None and not use_systemd_memory_limit:
        _warn_systemd_user_scope_unavailable()
        prlimit_args.append(f"--as={limits.mem_bytes}")
    if prlimit_args:
        cmd_list = ["prlimit", *prlimit_args, "--", *cmd_list]

    if use_systemd_memory_limit:
        cmd_list = [
            "systemd-run",
            "--user",
            "--scope",
            "--quiet",
            "-p",
            f"MemoryMax={limits.mem_bytes}",
            "-p",
            "MemorySwapMax=0",
            *cmd_list,
        ]

    return cmd_list


def _cleanup_process(proc: subprocess.Popen[bytes]) -> None:
    for pipe in (proc.stdin, proc.stdout, proc.stderr):
        if pipe is not None:
            with suppress(OSError, ValueError):
                pipe.close()
    if proc.poll() is None:
        with suppress(ProcessLookupError):
            proc.kill()
    proc.wait()


def _parse_gnu_time(
    timing_file: IO[bytes], ret_code: int
) -> tuple[float, float, float, float, int]:
    timing_file.seek(0)
    timing_lines = timing_file.read().decode("utf-8").strip().splitlines()
    timing_fields = timing_lines[-1].split() if timing_lines else []
    if len(timing_fields) != 5:
        raise RuntimeError(f"could not parse GNU time output: {timing_lines!r}")
    try:
        real_sec, user_sec, sys_sec = (float(value) for value in timing_fields[:3])
        maxrss_kb = float(timing_fields[3])
        command_ret_code = int(timing_fields[4])
    except ValueError as exc:
        raise RuntimeError(f"could not parse GNU time output: {timing_lines!r}") from exc

    signal_number = ret_code - 128
    if command_ret_code == 0 and 0 < signal_number < signal.NSIG:
        ret_code = -signal_number
    return real_sec, user_sec, sys_sec, maxrss_kb, ret_code


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

    # Resource limits deliberately wrap only the benchmarked command. Runner setup and
    # output compression stay outside them so benchmark limits do not depend on harness overhead.
    cmd_list = _cmd_list_with_limits(cmd, limits)

    env = os.environ.copy()
    if extra_env is not None:
        env.update(extra_env)

    with ExitStack() as stack:
        timing_file = stack.enter_context(tempfile.NamedTemporaryFile(mode="w+b"))
        command_stderr_file: IO[bytes] | None = None
        if limits.time_sec is not None:
            # Keep the command's stderr separate so timeout --verbose is an unambiguous
            # side channel. The shell immediately execs the command and leaves no supervisor.
            command_stderr_file = stack.enter_context(tempfile.NamedTemporaryFile(mode="w+b"))
            cmd_list = [
                "timeout",
                "--foreground",  # Let the command receive Ctrl-C from the terminal.
                "--verbose",
                "--kill-after=1",
                str(limits.time_sec),
                "/bin/sh",
                "-c",
                'exec 2>"$1"; shift; exec "$@"',
                "sh",
                command_stderr_file.name,
                *cmd_list,
            ]
        timed_cmd = (
            "/usr/bin/time",
            "-f",
            "%e %U %S %M %x",
            "-o",
            timing_file.name,
            "--",
            *cmd_list,
        )
        stdout: Any
        zstd_proc: subprocess.Popen[bytes] | None = None
        stdout_file: Any = None
        if stdout_path is not None:
            if compress_stdout:
                stdout_file = stack.enter_context(stdout_path.open("wb"))
                zstd_proc = subprocess.Popen(
                    ["zstd", "--fast", "-q"],  # noqa: S607
                    stdin=subprocess.PIPE,
                    stdout=stdout_file,
                )
                stack.callback(_cleanup_process, zstd_proc)
                stdout = zstd_proc.stdin
            else:
                stdout = stack.enter_context(stdout_path.open("wb"))
        else:
            stdout = subprocess.PIPE if stdout_to_str else subprocess.DEVNULL

        CMD_STR_LIM = 500
        cmd_str = " ".join(cmd)
        if len(cmd_str) > CMD_STR_LIM:
            cmd_str = cmd_str[: CMD_STR_LIM // 2] + "..." + cmd_str[-CMD_STR_LIM // 2 :]
        print(
            f"try_cmd: {cmd_str}, cwd: {cwd}, extra_env: {extra_env}, "
            f"stdout_to_str: {stdout_to_str}, stdout_path: {stdout_path}, "
            f"compress_stdout: {compress_stdout}, limits: {limits}"
        )
        try:
            completed = subprocess.run(
                timed_cmd,
                input=stdin_inp,
                stdout=stdout,
                stderr=subprocess.PIPE,
                cwd=cwd,
                env=env,
                check=False,
            )
        finally:
            # Close parent's copy of zstd's stdin so zstd sees EOF when the command exits.
            if zstd_proc is not None and zstd_proc.stdin is not None:
                zstd_proc.stdin.close()

        if zstd_proc is not None:
            zstd_proc.wait()
            assert stdout_file is not None
            stdout_file.close()

        stdout_str = completed.stdout.decode("utf-8") if completed.stdout else ""
        runner_stderr = completed.stderr or b""
        # GNU timeout returns 124 after TERM, or 137 if --kill-after needs SIGKILL.
        timed_out = (
            command_stderr_file is not None
            and completed.returncode in (124, 137)
            and bool(runner_stderr)
        )
        if command_stderr_file is None:
            stderr_bytes = runner_stderr
        else:
            command_stderr_file.seek(0)
            stderr_bytes = command_stderr_file.read()
            if not timed_out:
                stderr_bytes += runner_stderr
        stderr_str = stderr_bytes.decode("utf-8")

        if timed_out:
            raise TimeoutError(
                f"Command exceeded wall-clock time limit of {limits.time_sec} seconds: {cmd}"
            )

        real_sec, user_sec, sys_sec, maxrss_kb, ret_code = _parse_gnu_time(
            timing_file, completed.returncode
        )

        if zstd_proc is not None and zstd_proc.returncode != 0:
            raise RuntimeError(f"zstd compression failed with return code {zstd_proc.returncode}")

        if stdout_path is not None:
            if zstd_proc is None:
                assert not isinstance(stdout, (int, type(None)))
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
    env = os.environ.copy()
    if extra_env is not None:
        env.update(extra_env)
    print(f"cmd: {cmd}, cwd: {cwd}")
    ret_code = subprocess.run(["/bin/bash", "-c", cmd], cwd=cwd, env=env, check=False).returncode
    if ret_code != 0:
        click.echo(f"Running `{cmd}' failed with ret code {ret_code}.")
        raise RuntimeError(f"Shell command failed: {cmd}")
