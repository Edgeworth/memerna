# Copyright 2022 Eliot Courtney.
import copy
from dataclasses import dataclass
from dataclasses import field
from enum import Enum
from itertools import cycle
from itertools import islice
from pathlib import Path
import shutil

import click
from scripts.build.config import BuildCfg
from scripts.build.config import Sanitizer

AFL_MEMORY_LIMIT_MB = "2000"
AFL_TIME_LIMIT_MS = "2000"


class FuzzKind(str, Enum):
    REGULAR = "regular"
    ASAN = "asan"
    UBSAN = "ubsan"
    TSAN = "tsan"
    CFISAN = "cfisan"
    LAF = "laf"
    CMPLOG = "cmplog"

    def env(self) -> dict[str, str]:
        if self == FuzzKind.ASAN:
            return {"AFL_USE_ASAN": "1"}
        if self == FuzzKind.UBSAN:
            return {"AFL_USE_UBSAN": "1"}
        if self == FuzzKind.TSAN:
            return {"AFL_USE_TSAN": "1"}
        if self == FuzzKind.CFISAN:
            return {"AFL_USE_CFISAN": "1"}
        if self == FuzzKind.LAF:
            return {"AFL_LLVM_LAF_ALL": "1"}
        if self == FuzzKind.CMPLOG:
            return {"AFL_LLVM_CMPLOG": "1"}
        return {"AFL_HARDEN": "1 "}


@dataclass
class FuzzCfg:
    build_cfg: BuildCfg
    kind: FuzzKind = FuzzKind.REGULAR
    # extra args for afl-fuzz. not included in ident
    afl_args: list[str] = field(default_factory=list)
    index: int = 0  # which fuzzer this is when running multiple fuzzers

    def __post_init__(self) -> None:
        if self.error():
            raise ValueError(self.error())
        # Copy config so we can add our environment variables to it.
        self.build_cfg = copy.deepcopy(self.build_cfg)
        self.build_cfg.env.update(self.kind.env())

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

    def _build_single(self) -> None:
        # Build with this fuzz configs env vars.
        click.echo(f"Building fuzz configuration {self.ident()}")
        self.build_cfg.build(["fuzz"])

        # Copy artifacts to fuzz directory.
        self.bin_path().mkdir(parents=True, exist_ok=True)
        shutil.copy(self.build_cfg.build_path() / "fuzz", self.bin_path())

    # Note that this can't be called in parallel.
    def build(self) -> None:
        self._build_single()
        # If it was CMPLOG, put existing artifact as .cmplog.
        if self.kind == FuzzKind.CMPLOG:
            shutil.copy(self.bin_path() / "fuzz", self.bin_path() / "fuzz.cmplog")

    def _fuzz_cmd(self) -> str:
        cmd = f"./fuzz -md {self.build_cfg.src}/data "
        if self.build_cfg.rnastructure:
            cmd += f"-rd {self.build_cfg.src}/extern/miles_rnastructure/data_tables/ "
            cmd += "--mfe-rnastructure "  # "--subopt-rnastructure --part-rnastructure "

        cmd += "--afl"
        return cmd

    def afl_fuzz_cmd(self) -> str:
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
        cmd += " ".join(self.afl_args) + " "
        cmd += "-- " + self._fuzz_cmd()

        return cmd

    def afl_tmin_cmd(self, path: Path) -> str:
        cmd = ""

        cmd += f"afl-tmin -m {AFL_MEMORY_LIMIT_MB} -t {AFL_TIME_LIMIT_MS} "
        cmd += f"-i {path} -o {self.data_path()}/{path.name}.min "
        cmd += "-- " + self._fuzz_cmd()

        return cmd


def build_fuzz_cfgs(build_cfg: BuildCfg, max_num: int) -> list[FuzzCfg]:
    """Build an ensemble of fuzz configurations for a build configuration."""
    cfgs = []
    # Constructs fuzzers in this order:
    # 0: Regular fuzzer.
    # 1: LAF fuzzer.
    # 2: LAF fuzzer.
    # 3: LAF fuzzer.
    # 4: CMPLOG fuzzer.
    # 5: CMPLOG fuzzer with -l AT transformations.
    # 6: ASAN fuzzer.
    # 7: UBSAN fuzzer.
    # 8: TSAN fuzzer.
    # 9: CFISAN fuzzer.
    # 10: Regular fuzzers with combinations of -L 0, -Z, and different power schedules.
    kinds_args: list[tuple[FuzzKind, list[str]]] = [(FuzzKind.REGULAR, [])]

    # Skip other configurations if RNAstructure is enabled because compiling takes
    # too long.
    if not build_cfg.rnastructure:
        kinds_args += [
            (FuzzKind.CMPLOG, []),
            (FuzzKind.LAF, []),
            (FuzzKind.CMPLOG, ["-l", "AT"]),
            (FuzzKind.LAF, []),
            (FuzzKind.LAF, []),
            (FuzzKind.ASAN, []),
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
        cfgs.append(FuzzCfg(build_cfg=build_cfg, kind=kind, afl_args=extra_args, index=i))

    return cfgs[:max_num]
