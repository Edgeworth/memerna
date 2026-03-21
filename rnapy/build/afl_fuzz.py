# Copyright 2022 Eliot Courtney.
import copy
import dataclasses
import os
import shlex
import shutil
import subprocess
from dataclasses import dataclass, field
from enum import StrEnum
from functools import cached_property
from itertools import cycle, islice
from pathlib import Path

import click

from rnapy.build.build_cfg import BuildCfg, Sanitizer
from rnapy.model.model_cfg import CtdCfg, LonelyPairs

AFL_MEMORY_LIMIT_MB = "10000"
AFL_TIME_LIMIT_MS = "5000"
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
            return {"AFL_USE_CFISAN": "1", "AFL_CFISAN_VERBOSE": "1"}
        if self == AflFuzzKind.LAF:
            return {"AFL_LLVM_LAF_ALL": "1"}
        if self == AflFuzzKind.CMPLOG:
            return {"AFL_LLVM_CMPLOG": "1"}
        return {"AFL_HARDEN": "1"}


@dataclass
class AflFuzzCfg:
    build_cfg: BuildCfg

    fuzz_max_len: int | None = None
    fuzz_random_pseudofree: bool | None = None
    fuzz_energy_model: str | None = None
    fuzz_ctd: CtdCfg | None = None
    fuzz_lonely_pairs: LonelyPairs | None = None
    fuzz_backends: list[str] | None = None
    fuzz_brute_max: int | None = None
    fuzz_mfe: bool | None = None
    fuzz_mfe_rnastructure: bool | None = None
    fuzz_mfe_table: bool | None = None
    fuzz_subopt: bool | None = None
    fuzz_subopt_rnastructure: bool | None = None
    fuzz_subopt_strucs: int | None = None
    fuzz_subopt_delta: float | None = None
    fuzz_pfn: bool | None = None
    fuzz_pfn_rnastructure: bool | None = None

    kind: AflFuzzKind = AflFuzzKind.REGULAR
    # extra args for afl-fuzz. not included in ident
    afl_args: list[str] = field(default_factory=list)
    disable_trim: bool = False
    index: int = 0  # which fuzzer this is when running multiple fuzzers

    def __post_init__(self) -> None:
        if self.error():
            raise ValueError(self.error())
        # Keep the original build config to distinguish a set of fuzzers from another set.
        self.base_build_cfg = copy.deepcopy(self.build_cfg)
        # Copy config so we can add our environment variables to it.
        self.build_cfg = copy.deepcopy(self.build_cfg)
        # CFISAN trips in logger startup code, so force logging off only for that build.
        if self.kind == AflFuzzKind.CFISAN:
            self.build_cfg.enable_logging = False
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

    def _copy_built_binary(self, build_cfg: BuildCfg, dst_name: str) -> None:
        self.bin_path().mkdir(parents=True, exist_ok=True)
        shutil.copy(build_cfg.build_path() / AFL_TARGET, self.bin_path() / dst_name)

    # Note that this can't be called in parallel.
    def build(self) -> None:
        if self.kind != AflFuzzKind.CMPLOG:
            self._build_single()
            return

        click.echo(f"Building fuzz configuration {self.ident()} (regular + cmplog)")
        self.base_build_cfg.build([AFL_TARGET])
        self._copy_built_binary(self.base_build_cfg, AFL_TARGET)

        self.build_cfg.build([AFL_TARGET])
        self._copy_built_binary(self.build_cfg, AFL_TARGET + ".cmplog")

    @cached_property
    def afl_map_size(self) -> int:
        env = os.environ.copy()
        env["AFL_DUMP_MAP_SIZE"] = "1"
        res = subprocess.run(
            [f"./{AFL_TARGET}"],
            cwd=self.bin_path(),
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
        if res.returncode != 255:
            raise RuntimeError(
                f"Failed to query AFL_MAP_SIZE from {self.bin_path() / AFL_TARGET}: {res.stderr}"
            )
        return int(res.stdout.strip())

    def _afl_env(self) -> str:
        env = (
            "AFL_AUTORESUME=1 AFL_IMPORT_FIRST=1 AFL_TESTCACHE_SIZE=500 "
            "AFL_SKIP_CPUFREQ=1 AFL_CMPLOG_ONLY_NEW=1 "
        )
        env += f"AFL_MAP_SIZE={self.afl_map_size} "
        if self.disable_trim:
            env += "AFL_DISABLE_TRIM=1 "
        # Use ASAN options to enforce memory limit.
        # These use the AFL default ASAN options plus hard_rss_limit_mb
        if self.kind == AflFuzzKind.ASAN:
            env += (
                "ASAN_OPTIONS=abort_on_error=1:detect_leaks=0:malloc_context_size=0:"
                f"symbolize=0:allocator_may_return_null=1:hard_rss_limit_mb={AFL_MEMORY_LIMIT_MB} "
            )
        if self.kind == AflFuzzKind.TSAN:
            env += f"TSAN_OPTIONS=hard_rss_limit_mb={AFL_MEMORY_LIMIT_MB} "
        return env

    def _afl_limits(self) -> str:
        # ASAN allocates virtual memory which doesn't work well with AFL memory limit.
        if self.kind in [AflFuzzKind.ASAN, AflFuzzKind.TSAN]:
            return f"-t {AFL_TIME_LIMIT_MS}"
        return f"-m {AFL_MEMORY_LIMIT_MB} -t {AFL_TIME_LIMIT_MS}"

    def fuzz_argv(self) -> list[str]:
        cmd = [f"./{AFL_TARGET}"]
        if self.build_cfg.rnastructure:
            cmd += ["-rd", str(self.build_cfg.src / "extern/rnastructure_bridge/data_tables")]
        cmd += ["--memerna-data", str(self.build_cfg.src / "data")]
        if self.fuzz_max_len is not None:
            cmd += ["--max-len", str(self.fuzz_max_len)]
        if self.fuzz_random_pseudofree:
            cmd += ["--random-pf"]
        if self.fuzz_energy_model is not None:
            cmd += ["--energy-model", self.fuzz_energy_model]
        if self.fuzz_ctd is not None:
            cmd += ["--ctd", str(self.fuzz_ctd)]
        if self.fuzz_lonely_pairs is not None:
            cmd += ["--lonely-pairs", str(self.fuzz_lonely_pairs)]
        if self.fuzz_backends is not None:
            cmd += ["--backends", ",".join(self.fuzz_backends)]
        if self.fuzz_brute_max is not None:
            cmd += ["--brute-max", str(self.fuzz_brute_max)]
        if self.fuzz_mfe is not None:
            cmd += ["--mfe" if self.fuzz_mfe else "--no-mfe"]
        if self.fuzz_mfe_rnastructure is not None:
            cmd += ["--mfe-rnastructure" if self.fuzz_mfe_rnastructure else "--no-mfe-rnastructure"]
        if self.fuzz_mfe_table is not None:
            cmd += ["--mfe-table" if self.fuzz_mfe_table else "--no-mfe-table"]
        if self.fuzz_subopt is not None:
            cmd += ["--subopt" if self.fuzz_subopt else "--no-subopt"]
        if self.fuzz_subopt_rnastructure is not None:
            cmd += [
                "--subopt-rnastructure"
                if self.fuzz_subopt_rnastructure
                else "--no-subopt-rnastructure"
            ]
        if self.fuzz_subopt_strucs is not None:
            cmd += ["--subopt-strucs", str(self.fuzz_subopt_strucs)]
        if self.fuzz_subopt_delta is not None:
            cmd += ["--subopt-delta", str(self.fuzz_subopt_delta)]
        if self.fuzz_pfn is not None:
            cmd += ["--pfn" if self.fuzz_pfn else "--no-pfn"]
        if self.fuzz_pfn_rnastructure is not None:
            cmd += ["--pfn-rnastructure" if self.fuzz_pfn_rnastructure else "--no-pfn-rnastructure"]
        return cmd

    def fuzz_cmd(self) -> str:
        return " ".join(shlex.quote(arg) for arg in self.fuzz_argv())

    def afl_fuzz_cmd(self) -> str:
        cmd = ""

        instance = f"-M {self.ident()}" if self.index == 0 else f"-S {self.ident()}"

        # Add environment vars.
        cmd += f"{self._afl_env()} "
        # Add dictionary for fuzzing.
        afl_data_dir = self.build_cfg.src / AFL_DATA
        cmd += f"afl-fuzz -x {afl_data_dir}/dict.dct -i {afl_data_dir}/testcases "
        cmd += f"{self._afl_limits()} "
        cmd += f"-o {self.data_path()}/afl {instance} "
        if self.kind == AflFuzzKind.CMPLOG:
            cmd += f"-c ./{AFL_TARGET}.cmplog "
        cmd += " ".join(self.afl_args) + " "
        cmd += "-- " + self.fuzz_cmd()

        return cmd

    def afl_tmin_cmd(self, path: Path) -> str:
        cmd = ""

        cmd += f"AFL_MAP_SIZE={self.afl_map_size} "
        cmd += f"afl-tmin {self._afl_limits()}  "
        cmd += f"-i {path} -o {self.data_path()}/{path.name}.min "
        cmd += "-- " + self.fuzz_cmd()

        return cmd


def afl_fuzz_cfgs(afl_cfg: AflFuzzCfg, max_num_procs: int) -> list[AflFuzzCfg]:
    """Build an ensemble of fuzz configurations for a build configuration."""
    cfgs = []
    # Ensemble following AFL++ best practices (docs/fuzzing_in_depth.md):
    # 1. Main fuzzer (-M, gets old queue selection and no trimming automatically).
    # 2. CMPLOG with -l 2 (standard comparison logging).
    # 3. CMPLOG with -l 2AT (arithmetic + transformational solving).
    # 4. Sanitizer fuzzers (ASAN, UBSAN, TSAN, CFISAN).
    # 5. Extra regular fuzzers with varied power schedules, ~10% MOpt,
    #    ~10% old queue cycling (-Z), ~50% AFL_DISABLE_TRIM.
    # Note: LAF is disabled since it seems to cause heisenbugs.
    kinds_args: list[tuple[AflFuzzKind, list[str], bool]] = [
        (AflFuzzKind.REGULAR, [], False),
        (AflFuzzKind.CMPLOG, ["-l", "2"], False),
        (AflFuzzKind.CMPLOG, ["-l", "2AT"], False),
    ]

    # Skip sanitizer configurations if RNAstructure is enabled because we don't
    # care about these kinds of issues in RNAstructure.
    if not afl_cfg.build_cfg.rnastructure:
        kinds_args += [
            (AflFuzzKind.ASAN, [], False),
            (AflFuzzKind.UBSAN, [], False),
            (AflFuzzKind.TSAN, [], False),
            (AflFuzzKind.CFISAN, [], False),
        ]

    # ~10% MOpt (not well maintained but still recommended for a minority).
    mopt_cycle = [["-L", "0"]] + [[]] * 9
    # ~10% old queue cycling.
    queue_cycle = [["-Z"]] + [[]] * 9
    # ~50% disable trimming (recommended by AFL++ docs).
    trim_cycle = [True, False]
    # Power schedules: majority fast/explore, rest spread across recommended set.
    # Recommended by AFL++ docs: explore, fast, coe, lin, quad, exploit, rare.
    power_cycle = [
        ["-p", "fast"],
        ["-p", "fast"],
        ["-p", "fast"],
        ["-p", "explore"],
        ["-p", "explore"],
        ["-p", "exploit"],
        ["-p", "coe"],
        ["-p", "lin"],
        ["-p", "quad"],
        ["-p", "rare"],
    ]

    gen = zip(cycle(mopt_cycle), cycle(queue_cycle), cycle(power_cycle), cycle(trim_cycle))
    for mopt, queue, power, trim in islice(gen, max_num_procs):
        kinds_args.append((AflFuzzKind.REGULAR, mopt + queue + power, trim))

    for i, (kind, extra_args, trim) in enumerate(kinds_args):
        cfg_kwargs = {
            field.name: copy.deepcopy(getattr(afl_cfg, field.name))
            for field in dataclasses.fields(AflFuzzCfg)
        }
        cfg_kwargs["build_cfg"] = copy.deepcopy(afl_cfg.base_build_cfg)
        cfg_kwargs["kind"] = kind
        cfg_kwargs["afl_args"] = extra_args
        cfg_kwargs["disable_trim"] = trim
        cfg_kwargs["index"] = i
        cfg = AflFuzzCfg(**cfg_kwargs)
        cfgs.append(cfg)

    return cfgs[:max_num_procs]


def afl_fuzz_cfg_by_index(afl_cfg: AflFuzzCfg, index: int) -> AflFuzzCfg:
    if index < 0:
        raise ValueError("Fuzzer index must be >= 0.")
    return afl_fuzz_cfgs(afl_cfg, index + 1)[index]
