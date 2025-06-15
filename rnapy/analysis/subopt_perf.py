# Copyright 2022 Eliot Courtney.
import copy
from decimal import Decimal
from pathlib import Path

import click
import pandas as pd

from rnapy.bridge.memerna import MemeRna
from rnapy.bridge.rnapackage import RnaPackage
from rnapy.bridge.rnastructure import RNAstructure
from rnapy.bridge.viennarna import ViennaRna
from rnapy.data.memevault import MemeVault
from rnapy.model.model_cfg import CtdCfg, EnergyCfg, LonelyPairs, SuboptCfg
from rnapy.model.rna import Rna
from rnapy.util.util import row_by_key, strict_merge


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
    def _deltas() -> list[Decimal]:
        deltas = (
            "0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.5,1.8,2,3,4,5,6,7,8,9,10,"
            "15,20,25,30,35,45,60,75,100"
        )
        return sorted({Decimal(d) for d in deltas.split(",")})

    @staticmethod
    def _num_strucs() -> list[int]:
        return [1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000]

    def _run_once(
        self, program: RnaPackage, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg, rna: Rna
    ) -> bool:
        desc = program.desc(energy_cfg=energy_cfg, subopt_cfg=subopt_cfg)
        click.echo(f"Benchmarking folding with {desc} on {self.memevault.dataset}")
        output_path = self.output_dir / "subopt.json"

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

            row = row_by_key(output_path, data_keys)
            if row is not None:
                if row["failed"]:
                    click.echo(f"Skipping run {row} as it failed previously.")
                    return False
                click.echo(f"Skipping run {row} as it already exists in {output_path}")
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
            df = pd.DataFrame(data, index=[0])
            df.to_json(output_path, orient="records", lines=True, mode="a")

            if failed:
                return False
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
