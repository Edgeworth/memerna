# Copyright 2022 Eliot Courtney.
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
    delta: bool
    programs: list[tuple[RnaPackage, EnergyCfg, SuboptCfg, str]]

    def __init__(
        self,
        *,
        num_tries: int,
        time_sec_limit: int | None,
        mem_bytes_limit: int | None,
        memevault: MemeVault,
        output_dir: Path,
        delta: bool,
        memerna: MemeRna,
        rnastructure: RNAstructure,
        viennarna: ViennaRna,
    ) -> None:
        self.num_tries = num_tries
        self.memevault = memevault
        self.output_dir = output_dir
        self.delta = delta
        self.programs = [
            (
                rnastructure,
                EnergyCfg(ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.HEURISTIC),
                SuboptCfg(sorted_strucs=True, count_only=True),
                rnastructure.name(),
            ),
            (
                viennarna,
                EnergyCfg(ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.HEURISTIC),
                SuboptCfg(sorted_strucs=True, count_only=True),
                viennarna.name() + "-d3",
            ),
            (
                viennarna,
                EnergyCfg(ctd=CtdCfg.D2, lonely_pairs=LonelyPairs.HEURISTIC),
                SuboptCfg(sorted_strucs=True, count_only=True),
                viennarna.name() + "-d2",
            ),
            (
                memerna,
                EnergyCfg(ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.HEURISTIC),
                SuboptCfg(sorted_strucs=True, count_only=True, algorithm="iterative"),
                memerna.name() + "-delta-iterative",
            ),
            (
                memerna,
                EnergyCfg(ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.HEURISTIC),
                SuboptCfg(sorted_strucs=True, count_only=True, algorithm="persistent"),
                memerna.name() + "-delta-persistent",
            ),
        ]
        for program, _, _, _ in self.programs:
            program.limits.mem_bytes = mem_bytes_limit
            program.limits.time_sec = time_sec_limit

    @staticmethod
    def _max_delta(length: int) -> Decimal:
        max_delta = 3000 * Decimal(length) ** Decimal("-1.06")
        return SuboptPerfRunner._trunc_delta(max_delta)

    @staticmethod
    def _trunc_delta(delta: Decimal) -> Decimal:
        # Round down to nearest 0.1
        return Decimal(int(delta * 10)) / Decimal(10)

    @staticmethod
    def _deltas(length: int, count: int) -> list[Decimal]:
        max_delta = SuboptPerfRunner._max_delta(length)
        space = np.linspace(0.1, float(max_delta), count)
        values = [SuboptPerfRunner._trunc_delta(value) for value in space]
        return sorted(set(values))

    def _run_once(
        self,
        program: RnaPackage,
        energy_cfg: EnergyCfg,
        subopt_cfg: SuboptCfg,
        output_path: Path,
        rna: Rna,
    ) -> bool:
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
                    "delta": subopt_cfg.delta,
                    "num": subopt_cfg.strucs,
                    "num_strucs": rna_count,
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

    def run(self) -> None:
        for program, energy_cfg, subopt_cfg, name in self.programs:
            dataset = self.memevault.dataset
            click.echo(f"Benchmarking folding with {name} on {dataset}")
            kind = "delta" if self.delta else "num-strucs"
            output_path = self.output_dir / f"{dataset}_{name}_{kind}.results"
            if output_path.exists():
                click.echo(f"Output path {output_path} already exists, skipping")
                continue

            for rna_idx, rna in enumerate(self.memevault):
                deltas = self._deltas(len(rna), count=5)
                for delta in deltas:
                    click.echo(f"Running {program} on {rna_idx} {rna.name}, delta {delta}")
                    subopt_cfg.delta = delta
                    # Stop looking at larger deltas once a program fails once.
                    if not self._run_once(program, energy_cfg, subopt_cfg, output_path, rna):
                        click.echo(f"Failed, skipping remaining runs at {rna.name} for {program}")
                        break
