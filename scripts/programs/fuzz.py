# Copyright 2022 Eliot Courtney.
from dataclasses import dataclass
from enum import Enum
from itertools import cycle
from itertools import islice
import multiprocessing
from pathlib import Path
import shutil
from typing import Any

import click
import cloup
from scripts.build.args import build_cfg_from_args
from scripts.build.args import build_cfg_options
from scripts.build.cmake import generate_cmake
from scripts.build.config import BuildCfg
from scripts.build.config import Sanitizer
from scripts.util.command import run_shell
from scripts.util.util import fn_args


AFL_MEMORY_LIMIT_MB = "2000"
AFL_TIME_LIMIT_MS = "2000"


class FuzzKind(str, Enum):
    REGULAR = "regular"
    ASAN = "asan"
    MSAN = "msan"
    UBSAN = "ubsan"
    TSAN = "tsan"
    CFISAN = "cfisan"
    LAF = "laf"
    CMPLOG = "cmplog"

    def env(self) -> str:
        if self == FuzzKind.ASAN:
            return "AFL_USE_ASAN=1"
        if self == FuzzKind.MSAN:
            return "AFL_USE_MSAN=1"
        if self == FuzzKind.UBSAN:
            return "AFL_USE_UBSAN=1"
        if self == FuzzKind.TSAN:
            return "AFL_USE_TSAN=1"
        if self == FuzzKind.CFISAN:
            return "AFL_USE_CFISAN=1"
        if self == FuzzKind.LAF:
            return "AFL_LLVM_LAF_ALL=1"
        if self == FuzzKind.CMPLOG:
            return "AFL_LLVM_CMPLOG=1"
        return "AFL_HARDEN=1 "


@dataclass
class FuzzCfg:
    build_cfg: BuildCfg
    kind: FuzzKind
    extra_args: list[str]  # extra args for afl-fuzz. not included in ident
    index: int

    def error(self) -> str:
        if not self.build_cfg.is_afl():
            return "Fuzzing only supported for AFL configurations."
        if self.build_cfg.sanitizer != Sanitizer.NONE:
            return "Fuzzing not supported with build specified sanitizers."
        return ""

    def ident(self) -> str:
        return f"{self.index}-{self.kind}"

    def data_path(self) -> Path:
        """Directory where fuzzing data is stored."""
        return self.build_cfg.prefix / "memerna-fuzz" / self.build_cfg.ident()

    def bin_path(self) -> Path:
        """Directory where fuzzing binaries are stored"""
        return self.data_path() / self.ident()

    def _build_single(self, kind: FuzzKind) -> None:
        # Build with this fuzz configs env vars.
        build_path = self.build_cfg.build_path()
        # Need to clean because we changed environment variables.
        run_shell("make clean", cwd=build_path)

        build_cmd = f"{kind.env()} make -j$(($(nproc)-1)) fuzz"
        click.echo(f"Building fuzz configuration {build_cmd}")
        run_shell(build_cmd, cwd=build_path)

        # Copy artifacts to fuzz directory.
        self.bin_path().mkdir(parents=True, exist_ok=True)
        shutil.copy(build_path / "fuzz", self.bin_path())

    # Note that this can't be called in parallel.
    def build(self) -> None:
        build_path = self.build_cfg.build_path()
        self._build_single(self.kind)

        # If it was CMPLOG, rebuild with regular and put existing artifact as .cmplog.
        if self.kind == FuzzKind.CMPLOG:
            shutil.copy(build_path / "fuzz", self.bin_path() / "fuzz.cmplog")
            self._build_single(FuzzKind.REGULAR)
            shutil.copy(build_path / "fuzz", self.bin_path())

    def fuzz_cmd(self) -> str:
        cmd = ""

        instance = f"-M {self.ident()}" if self.index == 0 else f"-S {self.ident()}"

        # Add environment vars.
        cmd += "AFL_AUTORESUME=1 AFL_IMPORT_FIRST=1 AFL_TESTCACHE_SIZE=500 AFL_SKIP_CPUFREQ=1 "
        # Add dictionary for fuzzing.
        cmd += f"afl-fuzz -x {self.build_cfg.src}/extern/aflplusplus/fuzz/dict.dct "
        cmd += f"-m {AFL_MEMORY_LIMIT_MB} -t {AFL_TIME_LIMIT_MS} "
        cmd += f"-i {self.build_cfg.src}/extern/aflplusplus/fuzz/testcases "
        cmd += f"-o {self.data_path()}/afl {instance} "
        if self.kind == FuzzKind.CMPLOG:
            cmd += "-c ./fuzz.cmplog "
        cmd += " ".join(self.extra_args) + " "
        cmd += "-- "

        cmd += f"./fuzz -md {self.build_cfg.src}/data "
        if self.build_cfg.rnastructure:
            cmd += f"-rd {self.build_cfg.src}/extern/miles_rnastructure/data_tables "
            cmd += "--mfe-rnastructure "  # "--subopt-rnastructure --part-rnastructure "

        cmd += "--afl "

        return cmd


def build_fuzz_cfgs(build_cfg: BuildCfg, max_num: int) -> list[FuzzCfg]:
    cfgs = []
    # Constructs fuzzers in this order:
    # 0: Regular fuzzer.
    # 1: LAF fuzzer.
    # 2: LAF fuzzer.
    # 3: LAF fuzzer.
    # 4: CMPLOG fuzzer.
    # 5: CMPLOG fuzzer with -l AT transformations.
    # 6: ASAN fuzzer.
    # 7: MSAN fuzzer.
    # 8: UBSAN fuzzer.
    # 9: TSAN fuzzer.
    # 10: CFISAN fuzzer.
    # 11: Regular fuzzers with combinations of -L 0, -Z, and different power schedules.
    kinds_args: list[tuple[FuzzKind, list[str]]] = [
        (FuzzKind.REGULAR, []),
        (FuzzKind.CMPLOG, []),
        (FuzzKind.LAF, []),
        (FuzzKind.CMPLOG, ["-l", "AT"]),
        (FuzzKind.LAF, []),
        (FuzzKind.LAF, []),
        (FuzzKind.ASAN, []),
        (FuzzKind.MSAN, []),
        (FuzzKind.UBSAN, []),
        (FuzzKind.TSAN, []),
        (FuzzKind.CFISAN, []),
    ]

    # 1/4th of the time.
    mopt_cycle = [["-L", "0"]] + [[]] * 3
    queue_cycle = [["-Z"]] + [[]] * 10
    power_cycle = [
        ["-p", "fast"],
        ["-p", "explore"],
        ["-p", "exploit"],
        ["-p", "seek"],
        ["-p", "rare"],
        ["-p", "mmopt"],
        ["-p", "coe"],
        ["-p", "lin"],
        ["-p", "quad"],
    ]

    gen = zip(cycle(mopt_cycle), cycle(queue_cycle), cycle(power_cycle))
    for mopt, queue, power in islice(gen, max_num):
        kinds_args.append((FuzzKind.REGULAR, mopt + queue + power))

    for i, (kind, extra_args) in enumerate(kinds_args):
        cfgs.append(FuzzCfg(build_cfg=build_cfg, kind=kind, extra_args=extra_args, index=i))

    return cfgs[:max_num]


def run_fuzz(cfg: FuzzCfg) -> None:
    cmd = cfg.fuzz_cmd()
    click.echo(f"Running fuzz {cmd}")
    run_shell(cmd, cwd=cfg.bin_path())


@cloup.command()
@build_cfg_options
@cloup.option(
    "--num-proc",
    default=multiprocessing.cpu_count() - 1,
    help="Number of fuzzing configurations to run.",
)
def fuzz(
    num_proc: int,
    **_kwargs: Any,
) -> None:
    build_cfg = build_cfg_from_args(**fn_args())
    cfgs = build_fuzz_cfgs(build_cfg, num_proc)

    for cfg in cfgs:
        if cfg.error():
            raise click.UsageError(cfg.error())

    generate_cmake(build_cfg)
    for cfg in cfgs:
        cfg.build()

    with multiprocessing.Pool(len(cfgs)) as pool:
        pool.map(run_fuzz, cfgs)
