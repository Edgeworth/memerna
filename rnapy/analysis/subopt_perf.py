# Copyright 2022 Eliot Courtney.
import copy
import threading
from decimal import Decimal
from pathlib import Path
from typing import Any

import click
import polars as pl

from rnapy.bridge.memerna import MemeRna
from rnapy.bridge.rnapackage import RnaPackage
from rnapy.bridge.rnastructure import RNAstructure
from rnapy.bridge.viennarna import ViennaRna
from rnapy.data.memevault import MemeVault
from rnapy.model.model_cfg import CtdCfg, EnergyCfg, LonelyPairs, SuboptCfg
from rnapy.model.rna import Rna
from rnapy.util.util import append_ndjson, parallel_map, strict_merge


class SuboptPerfRunner:
    num_tries: int
    memevault: MemeVault
    output_path: Path
    rna_lengths: tuple[int, ...]
    jobs: int
    rnas: list[Rna]
    programs: list[tuple[RnaPackage, EnergyCfg, SuboptCfg]]
    _file_lock: threading.Lock
    _cached: pl.DataFrame

    def __init__(
        self,
        *,
        num_tries: int,
        memevault: MemeVault,
        output_path: Path,
        rna_lengths: tuple[int, ...],
        jobs: int,
        memerna: MemeRna,
        rnastructure: RNAstructure,
        viennarna: ViennaRna,
    ) -> None:
        self.num_tries = num_tries
        self.memevault = memevault
        self.output_path = output_path
        self.rna_lengths = rna_lengths
        self.jobs = jobs
        self._file_lock = threading.Lock()
        # Pre-materialize RNAs since sqlite3 connections aren't thread-safe.
        self.rnas = [rna for rna in memevault if not rna_lengths or len(rna) in rna_lengths]
        # Load existing results once to avoid re-reading the file per run.
        if output_path.exists() and output_path.stat().st_size > 0:
            self._cached = pl.read_ndjson(output_path)
        else:
            self._cached = pl.DataFrame()
        self.programs = [
            (
                rnastructure,
                EnergyCfg(ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.HEURISTIC, energy_model="t04"),
                SuboptCfg(sorted_strucs=True, count_only=True),
            ),
            (
                viennarna,
                EnergyCfg(ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.HEURISTIC, energy_model="t04"),
                SuboptCfg(sorted_strucs=True, count_only=True),
            ),
            (
                viennarna,
                EnergyCfg(ctd=CtdCfg.D2, lonely_pairs=LonelyPairs.HEURISTIC, energy_model="t04"),
                SuboptCfg(sorted_strucs=True, count_only=True),
            ),
            (
                memerna,
                EnergyCfg(
                    ctd=CtdCfg.ALL,
                    lonely_pairs=LonelyPairs.HEURISTIC,
                    energy_model="t04",
                    backend="baseopt",
                ),
                SuboptCfg(sorted_strucs=True, count_only=True, algorithm="iterative"),
            ),
            (
                memerna,
                EnergyCfg(
                    ctd=CtdCfg.D2,
                    lonely_pairs=LonelyPairs.HEURISTIC,
                    energy_model="t04",
                    backend="base",
                ),
                SuboptCfg(sorted_strucs=True, count_only=True, algorithm="iterative"),
            ),
            (
                memerna,
                EnergyCfg(
                    ctd=CtdCfg.ALL,
                    lonely_pairs=LonelyPairs.HEURISTIC,
                    energy_model="t04",
                    backend="baseopt",
                ),
                SuboptCfg(sorted_strucs=True, count_only=True, algorithm="iterative-lowmem"),
            ),
            (
                memerna,
                EnergyCfg(
                    ctd=CtdCfg.D2,
                    lonely_pairs=LonelyPairs.HEURISTIC,
                    energy_model="t04",
                    backend="base",
                ),
                SuboptCfg(sorted_strucs=True, count_only=True, algorithm="iterative-lowmem"),
            ),
            (
                memerna,
                EnergyCfg(
                    ctd=CtdCfg.ALL,
                    lonely_pairs=LonelyPairs.HEURISTIC,
                    energy_model="t04",
                    backend="baseopt",
                ),
                SuboptCfg(sorted_strucs=True, count_only=True, algorithm="persistent"),
            ),
            (
                memerna,
                EnergyCfg(
                    ctd=CtdCfg.D2,
                    lonely_pairs=LonelyPairs.HEURISTIC,
                    energy_model="t04",
                    backend="base",
                ),
                SuboptCfg(sorted_strucs=True, count_only=True, algorithm="persistent"),
            ),
            (
                memerna,
                EnergyCfg(
                    ctd=CtdCfg.ALL,
                    lonely_pairs=LonelyPairs.HEURISTIC,
                    energy_model="t04",
                    backend="baseopt",
                ),
                SuboptCfg(sorted_strucs=True, count_only=True, algorithm="persistent-lowmem"),
            ),
            (
                memerna,
                EnergyCfg(
                    ctd=CtdCfg.D2,
                    lonely_pairs=LonelyPairs.HEURISTIC,
                    energy_model="t04",
                    backend="base",
                ),
                SuboptCfg(sorted_strucs=True, count_only=True, algorithm="persistent-lowmem"),
            ),
        ]

    @staticmethod
    def _deltas(rna_length: int) -> list[Decimal]:
        max_delta = 61 if rna_length < 500 else 31
        return [Decimal(i) / 10 for i in range(max_delta)]

    @staticmethod
    def _num_strucs() -> list[int]:
        return list(range(100_000, 10_000_001, 100_000))

    def _lookup_cached(self, data_keys: dict) -> dict[str, Any] | None:
        df = self._cached
        if df.is_empty():
            return None
        for key, value in data_keys.items():
            if key not in df.columns:
                return None
            if value is None:
                df = df.filter(pl.col(key).is_null())
            else:
                df = df.filter(pl.col(key) == value)
        if df.is_empty():
            return None
        return df.row(0, named=True)

    def _run_once(
        self, program: RnaPackage, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg, rna: Rna
    ) -> bool:
        desc = program.desc(energy_cfg=energy_cfg, subopt_cfg=subopt_cfg)
        click.echo(f"Benchmarking folding with {desc} on {self.memevault.dataset}")

        for run_idx in range(self.num_tries):
            data_keys = strict_merge(
                desc,
                {
                    "dataset": self.memevault.dataset,
                    "rna_name": rna.name,
                    "rna_length": len(rna),
                    "run_idx": run_idx,
                },
            )

            row = self._lookup_cached(data_keys)
            if row is not None:
                click.echo(f"Skipping run {data_keys} (cached)")
                if row["failed"]:
                    return False
                continue

            failed = False
            data_values: dict = {}
            rna_count = None
            try:
                rna_count, cmd_res = program.subopt(rna, energy_cfg, subopt_cfg)
                failed = cmd_res.ret_code != 0
                data_values = strict_merge(
                    data_values,
                    {
                        "output_strucs": rna_count,
                        "maxrss_bytes": cmd_res.maxrss_bytes,
                        "user_sec": cmd_res.user_sec,
                        "sys_sec": cmd_res.sys_sec,
                        "real_sec": cmd_res.real_sec,
                    },
                )
            except Exception as e:
                click.echo(f"Error running {program} on {rna.name}: {e}")
                failed = True

            assert failed or isinstance(rna_count, int), f"Expected int, got {type(rna_count)}"

            data = strict_merge(data_keys, data_values, {"failed": failed})
            with self._file_lock:
                append_ndjson(self.output_path, pl.DataFrame([data]))

            if failed:
                return False
        return True

    def _run_rna(
        self,
        program: RnaPackage,
        energy_cfg: EnergyCfg,
        subopt_cfgs: list[SuboptCfg],
        rna_idx: int,
        rna: Rna,
    ) -> None:
        for cfg in subopt_cfgs:
            click.echo(f"Running {program} on {rna_idx} {rna.name}, cfg {cfg}")
            if not self._run_once(program, energy_cfg, cfg, rna):
                click.echo(f"Failed, skipping remaining runs at {rna.name} for {program}")
                break

    def run(self) -> None:
        jobs: list[tuple[RnaPackage, EnergyCfg, list[SuboptCfg], int, Rna]] = []
        for program, energy_cfg, base_cfg in self.programs:
            strucs_cfgs = []
            for num_strucs in self._num_strucs():
                cfg = copy.deepcopy(base_cfg)
                cfg.strucs = num_strucs
                strucs_cfgs.append(cfg)

            for rna_idx, rna in enumerate(self.rnas):
                delta_cfgs = []
                for delta in self._deltas(len(rna)):
                    cfg = copy.deepcopy(base_cfg)
                    cfg.delta = delta
                    delta_cfgs.append(cfg)

                jobs.append((program, energy_cfg, delta_cfgs, rna_idx, rna))
                jobs.append((program, energy_cfg, strucs_cfgs, rna_idx, rna))

        parallel_map(self._run_rna, jobs, max_workers=self.jobs)
