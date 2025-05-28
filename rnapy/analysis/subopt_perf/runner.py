# Copyright 2022 Eliot Courtney.
import copy
from decimal import Decimal
from pathlib import Path

import click
import numpy as np
import pandas as pd

from rnapy.bridge.memerna import MemeRna
from rnapy.bridge.rnapackage import RnaPackage
from rnapy.bridge.rnastructure import RNAstructure
from rnapy.bridge.viennarna import ViennaRna
from rnapy.data.memevault import MemeVault
from rnapy.model.model_cfg import CtdCfg, EnergyCfg, LonelyPairs, SuboptCfg
from rnapy.model.rna import Rna


class SuboptPerfRunner:
    num_tries: int
    memevault: MemeVault
    output_dir: Path
    programs: list[tuple[RnaPackage, EnergyCfg, SuboptCfg]]

    def __init__(
        self,
        *,
        num_tries: int,
        memevault: MemeVault,
        output_dir: Path,
        memerna: MemeRna,
        rnastructure: RNAstructure,
        viennarna: ViennaRna,
    ) -> None:
        self.num_tries = num_tries
        self.memevault = memevault
        self.output_dir = output_dir
        self.programs = [
            (
                rnastructure,
                EnergyCfg(ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.HEURISTIC),
                SuboptCfg(sorted_strucs=True, count_only=True),
            ),
            (
                viennarna,
                EnergyCfg(ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.HEURISTIC),
                SuboptCfg(sorted_strucs=True, count_only=True),
            ),
            (
                viennarna,
                EnergyCfg(ctd=CtdCfg.D2, lonely_pairs=LonelyPairs.HEURISTIC),
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
    def _trunc_delta(delta: float) -> Decimal:
        # Round down to nearest 0.1
        return Decimal(int(delta * 10)) / Decimal(10)

    @staticmethod
    def _deltas() -> list[Decimal]:
        # We keep running until we time out.
        space = np.arange(0.1, 100.0, 0.1)
        values = [SuboptPerfRunner._trunc_delta(float(value)) for value in space]
        return sorted(set(values))

    @staticmethod
    def _num_strucs() -> list[int]:
        # We keep running until we time out.
        space = np.logspace(0, 9, num=10, base=10, dtype=int)
        return [int(value) for value in space]

    def _run_once(
        self, program: RnaPackage, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg, rna: Rna
    ) -> bool:
        desc = program.desc(energy_cfg=energy_cfg, subopt_cfg=subopt_cfg)
        dataset = self.memevault.dataset
        click.echo(f"Benchmarking folding with {desc} on {dataset}")
        output_path = self.output_dir / f"{dataset}_{desc}.results"
        if output_path.exists():
            click.echo(f"Output path {output_path} already exists, skipping")
            return True

        for run_idx in range(self.num_tries):
            try:
                rna_count, cmd_res = program.subopt(rna, energy_cfg, subopt_cfg)
            except Exception as e:
                click.echo(f"Error running {program} on {rna.name}: {e}")
                return False
            assert isinstance(rna_count, int), f"Expected int, got {type(rna_count)}"

            df = pd.DataFrame(
                {
                    "name": rna.name,
                    "run_idx": run_idx,
                    "length": len(rna),
                    "input_delta": subopt_cfg.delta,
                    "input_strucs": subopt_cfg.strucs,
                    "output_strucs": rna_count,
                    "maxrss_bytes": cmd_res.maxrss_bytes,
                    "user_sec": cmd_res.user_sec,
                    "sys_sec": cmd_res.sys_sec,
                    "real_sec": cmd_res.real_sec,
                },
                index=[0],
            )
            df.to_csv(
                output_path,
                mode="a",
                header=not output_path.exists(),
                index=False,
                float_format="%g",
            )
        return True

    def _run(
        self, program: RnaPackage, energy_cfg: EnergyCfg, subopt_cfgs: list[SuboptCfg]
    ) -> None:
        for rna_idx, rna in enumerate(self.memevault):
            for cfg in subopt_cfgs:
                click.echo(f"Running {program} on {rna_idx} {rna.name}, cfg {cfg}")
                # Stop looking at larger cfgs (e.g. larger deltas) once a program fails once.
                if not self._run_once(program, energy_cfg, cfg, rna):
                    click.echo(f"Failed, skipping remaining runs at {rna.name} for {program}")
                    break

    def run(self) -> None:
        for program, energy_cfg, base_cfg in self.programs:
            cfgs = []
            for delta in self._deltas():
                cfg = copy.deepcopy(base_cfg)
                cfg.delta = delta
                cfgs.append(cfg)

            self._run(program, energy_cfg, cfgs)

            cfgs = []
            for num_strucs in self._num_strucs():
                cfg = copy.deepcopy(base_cfg)
                cfg.strucs = num_strucs
                cfgs.append(cfg)

            self._run(program, energy_cfg, cfgs)
