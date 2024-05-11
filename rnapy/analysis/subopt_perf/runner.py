# Copyright 2022 Eliot Courtney.
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


class SuboptPerfRunner:
    num_tries: int
    memevault: MemeVault
    output_dir: Path
    delta: Decimal | None
    programs: list[tuple[RnaPackage, EnergyCfg, SuboptCfg, str]]

    def __init__(
        self,
        *,
        num_tries: int,
        time_sec_limit: int | None,
        mem_bytes_limit: int | None,
        memevault: MemeVault,
        output_dir: Path,
        delta: Decimal | None,
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
                SuboptCfg(delta=delta, sorted_strucs=True),
                rnastructure.name(),
            ),
            (
                viennarna,
                EnergyCfg(ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.HEURISTIC),
                SuboptCfg(delta=delta, sorted_strucs=True),
                viennarna.name() + "-d3",
            ),
            (
                viennarna,
                EnergyCfg(ctd=CtdCfg.D2, lonely_pairs=LonelyPairs.HEURISTIC),
                SuboptCfg(delta=delta, sorted_strucs=True),
                viennarna.name() + "-d2",
            ),
            (
                memerna,
                EnergyCfg(ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.HEURISTIC),
                SuboptCfg(delta=delta, sorted_strucs=True),
                memerna.name(),
            ),
        ]
        for program, _, _, _ in self.programs:
            program.limits.mem_bytes = mem_bytes_limit
            program.limits.time_sec = time_sec_limit

    def run(self) -> None:
        for program, energy_cfg, subopt_cfg, name in self.programs:
            dataset = self.memevault.dataset
            click.echo(f"Benchmarking folding with {name} on {dataset}")
            output_path = self.output_dir / f"{dataset}_{name}_{self.delta}.results"
            if output_path.exists():
                raise RuntimeError(f"Output path {output_path} already exists")

            for rna_idx, rna in enumerate(self.memevault):
                click.echo(f"Running {program} on {rna_idx} {rna.name}")
                failed = False
                for run_idx in range(self.num_tries):
                    try:
                        rnas, cmd_res = program.subopt(rna, energy_cfg, subopt_cfg)
                    except Exception as e:
                        click.echo(f"Error running {program} on {rna.name}: {e}")
                        failed = True
                        break
                    df = pd.DataFrame(
                        {
                            "name": rna.name,
                            "run_idx": run_idx,
                            "length": len(rnas),
                            "num_strucs": len(rnas),
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
                if failed:
                    click.echo(f"Failed, skipping remaining runs at {rna.name} for {program}")
                    break
                break
