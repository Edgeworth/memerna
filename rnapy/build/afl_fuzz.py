# Copyright 2022 Eliot Courtney.
import copy
import dataclasses
import shutil
from dataclasses import dataclass, field
from enum import StrEnum
from itertools import cycle, islice
from pathlib import Path

import click

from rnapy.build.build_cfg import BuildCfg, Sanitizer

AFL_MEMORY_LIMIT_MB = "2000"
AFL_TIME_LIMIT_MS = "2000"
AFL_TARGET = "fuzz_afl"
AFL_DATA = Path("data") / "aflplusplus" / AFL_TARGET


class AflFuzzKind(StrEnum):
    REGULAR = "regular"
    ASAN = "asan"
    UBSAN = "ubsan"
    TSAN = "tsan"
    CFISAN = "cfisan"
    LAF = "laf"
    CMPLOG = "cmplog"

    def env(self) -> dict[str, str]:
        if self == AflFuzzKind.ASAN:
            return {"AFL_USE_ASAN": "1"}
        if self == AflFuzzKind.UBSAN:
            return {"AFL_USE_UBSAN": "1"}
        if self == AflFuzzKind.TSAN:
            return {"AFL_USE_TSAN": "1"}
        if self == AflFuzzKind.CFISAN:
            return {"AFL_USE_CFISAN": "1"}
        if self == AflFuzzKind.LAF:
            return {"AFL_LLVM_LAF_ALL": "1"}
        if self == AflFuzzKind.CMPLOG:
            return {"AFL_LLVM_CMPLOG": "1"}
        return {"AFL_HARDEN": "1"}


@dataclass
class AflFuzzCfg:
    build_cfg: BuildCfg

    # Info for the actual fuzz_afl invocation:
    fuzz_max_len: int
    fuzz_seed: int | None  # For fixed random model
    fuzz_random_models: bool  # For random models every fuzz invocation
    fuzz_energy_models: list[str]
    fuzz_brute_max: int
    fuzz_mfe: bool
    fuzz_mfe_rnastructure: bool
    fuzz_mfe_table: bool
    fuzz_subopt: bool
    fuzz_subopt_rnastructure: bool
    fuzz_subopt_strucs: int
    fuzz_subopt_delta: float
    fuzz_part: bool
    fuzz_part_rnastructure: bool

    kind: AflFuzzKind = AflFuzzKind.REGULAR
    # extra args for afl-fuzz. not included in ident
    afl_args: list[str] = field(default_factory=list)
    index: int = 0  # which fuzzer this is when running multiple fuzzers

    def __post_init__(self) -> None:
        if self.error():
            raise ValueError(self.error())
        # Keep the original build config to distinguish a set of fuzzers from another set.
        self.base_build_cfg = copy.deepcopy(self.build_cfg)
        # Copy config so we can add our environment variables to it.
        self.build_cfg = copy.deepcopy(self.build_cfg)
        self.build_cfg.env.update(self.kind.env())

        # Remove mutually exclusive env vars.
        if self.kind == AflFuzzKind.ASAN:
            del self.build_cfg.env["AFL_HARDEN"]

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
        return self.build_cfg.prefix / "memerna-afl" / self.base_build_cfg.ident()

    def bin_path(self) -> Path:
        """Directory where fuzzing binaries are stored"""
        return self.data_path() / self.ident()

    def _build_single(self) -> None:
        # Build with this fuzz configs env vars.
        click.echo(f"Building fuzz configuration {self.ident()}")
        self.build_cfg.build([AFL_TARGET])

        # Copy artifacts to fuzz directory.
        self.bin_path().mkdir(parents=True, exist_ok=True)
        shutil.copy(self.build_cfg.build_path() / AFL_TARGET, self.bin_path())

    # Note that this can't be called in parallel.
    def build(self) -> None:
        self._build_single()
        # If it was CMPLOG, put existing artifact as .cmplog.
        if self.kind == AflFuzzKind.CMPLOG:
            shutil.copy(self.bin_path() / AFL_TARGET, self.bin_path() / (AFL_TARGET + ".cmplog"))

    def _afl_env(self) -> str:
        # Use ASAN options to enforce memory limit.
        # These use the AFL default ASAN options plus hard_rss_limit_mb
        if self.kind == AflFuzzKind.ASAN:
            return (
                "ASAN_OPTIONS=abort_on_error=1:detect_leaks=0:malloc_context_size=0:"
                f"symbolize=0:allocator_may_return_null=1:hard_rss_limit_mb={AFL_MEMORY_LIMIT_MB}"
            )
        if self.kind == AflFuzzKind.TSAN:
            return f"TSAN_OPTIONS=hard_rss_limit_mb={AFL_MEMORY_LIMIT_MB}"
        return ""

    def _afl_limits(self) -> str:
        # ASAN allocates virtual memory which doesn't work well with AFL memory limit.
        if self.kind in [AflFuzzKind.ASAN, AflFuzzKind.TSAN]:
            return f"-t {AFL_TIME_LIMIT_MS}"
        return f"-m {AFL_MEMORY_LIMIT_MB} -t {AFL_TIME_LIMIT_MS}"

    def _fuzz_cmd(self) -> str:
        cmd = f"./{AFL_TARGET} "
        if self.build_cfg.rnastructure:
            cmd += f"-rd {self.build_cfg.src}/extern/rnastructure_bridge/data_tables/ "
        cmd += f"--memerna-data {self.build_cfg.src / 'data'} "
        cmd += f"--max-len {self.fuzz_max_len} "
        if self.fuzz_seed:
            cmd += f"--seed {self.fuzz_seed} "
        if self.fuzz_random_models:
            cmd += "--random-models "
        cmd += f"--energy-models {','.join(self.fuzz_energy_models)} "
        cmd += f"--brute-max {self.fuzz_brute_max} "
        cmd += "--mfe " if self.fuzz_mfe else "--no-mfe "
        cmd += "--mfe-rnastructure " if self.fuzz_mfe_rnastructure else "--no-mfe-rnastructure "
        cmd += "--mfe-table " if self.fuzz_mfe_table else "--no-mfe-table "
        cmd += "--subopt " if self.fuzz_subopt else "--no-subopt "
        cmd += (
            "--subopt-rnastructure "
            if self.fuzz_subopt_rnastructure
            else "--no-subopt-rnastructure "
        )
        cmd += f"--subopt-strucs {self.fuzz_subopt_strucs} "
        cmd += f"--subopt-delta {self.fuzz_subopt_delta} "
        cmd += "--part " if self.fuzz_part else "--no-part-rnastructure "
        cmd += "--part-rnastructure " if self.fuzz_part_rnastructure else "--no-part-rnastructure "

        return cmd

    def afl_fuzz_cmd(self) -> str:
        cmd = ""

        instance = f"-M {self.ident()}" if self.index == 0 else f"-S {self.ident()}"

        # Add environment vars.
        cmd += "AFL_AUTORESUME=1 AFL_IMPORT_FIRST=1 AFL_TESTCACHE_SIZE=500 AFL_SKIP_CPUFREQ=1 "
        cmd += f"{self._afl_env()} "
        # Add dictionary for fuzzing.
        afl_data_dir = self.build_cfg.src / AFL_DATA
        cmd += f"afl-fuzz -x {afl_data_dir}/dict.dct -i {afl_data_dir}/testcases "
        cmd += f"{self._afl_limits()} "
        cmd += f"-o {self.data_path()}/afl {instance} "
        if self.kind == AflFuzzKind.CMPLOG:
            cmd += f"-c ./{AFL_TARGET}.cmplog "
        cmd += " ".join(self.afl_args) + " "
        cmd += "-- " + self._fuzz_cmd()

        return cmd

    def afl_tmin_cmd(self, path: Path) -> str:
        cmd = ""

        cmd += f"afl-tmin {self._afl_limits()}  "
        cmd += f"-i {path} -o {self.data_path()}/{path.name}.min "
        cmd += "-- " + self._fuzz_cmd()

        return cmd


def afl_fuzz_cfgs(afl_cfg: AflFuzzCfg, max_num_procs: int) -> list[AflFuzzCfg]:
    """Build an ensemble of fuzz configurations for a build configuration."""
    cfgs = []
    # Constructs fuzzers in this order:
    # Regular fuzzer - the main fuzzer.
    # CMPLOG fuzzer.
    # CMPLOG fuzzer with -l AT transformations.
    # ASAN fuzzer.
    # UBSAN fuzzer.
    # TSAN fuzzer.
    # CFISAN fuzzer.
    # Regular fuzzers with combinations of -L 0, -Z, and different power schedules.
    # Note: LAF is disabled since it seems to cause heisenbugs.
    kinds_args: list[tuple[AflFuzzKind, list[str]]] = [
        (AflFuzzKind.REGULAR, []),
        (AflFuzzKind.CMPLOG, []),
        (AflFuzzKind.CMPLOG, ["-l", "AT"]),
    ]

    # Skip other configurations if RNAstructure is enabled because we don't
    # care about these kinds of issues in RNAstructure.
    # Also skip LAF, as it takes too long to compile.
    if not afl_cfg.build_cfg.rnastructure:
        kinds_args += [
            (AflFuzzKind.ASAN, []),
            (AflFuzzKind.UBSAN, []),
            (AflFuzzKind.TSAN, []),
            (AflFuzzKind.CFISAN, []),
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
    for mopt, queue, power in islice(gen, max_num_procs):
        kinds_args.append((AflFuzzKind.REGULAR, mopt + queue + power))

    for i, (kind, extra_args) in enumerate(kinds_args):
        cfg = dataclasses.replace(afl_cfg, kind=kind, afl_args=extra_args, index=i)
        cfgs.append(cfg)

    return cfgs[:max_num_procs]
