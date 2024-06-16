# Copyright 2022 Eliot Courtney.
from pathlib import Path

import click
import pandas as pd

from rnapy.bridge.linearfold import LinearFold
from rnapy.bridge.memerna01 import MemeRna01
from rnapy.bridge.rnapackage import RnaPackage
from rnapy.bridge.rnastructure import RNAstructure
from rnapy.bridge.sparsemfefold import SparseMFEFold
from rnapy.bridge.sparsernafold import SparseRNAFolD
from rnapy.bridge.viennarna import ViennaRna
from rnapy.data.memevault import MemeVault
from rnapy.model.model_cfg import CtdCfg, EnergyCfg, LonelyPairs


class FoldPerfRunner:
    num_tries: int
    memevault: MemeVault
    output_dir: Path
    programs: list[tuple[RnaPackage, EnergyCfg]]

    def __init__(
        self,
        *,
        num_tries: int,
        time_sec_limit: int | None,
        mem_bytes_limit: int | None,
        memevault: MemeVault,
        output_dir: Path,
        memerna01: MemeRna01,
        linearfold: LinearFold,
        rnastructure: RNAstructure,
        viennarna: ViennaRna,
        sparsemfefold: SparseMFEFold,
        sparsernafold: SparseRNAFolD,
    ) -> None:
        self.num_tries = num_tries
        self.memevault = memevault
        self.output_dir = output_dir
        self.programs = [
            (
                memerna01,
                EnergyCfg(
                    ctd=CtdCfg.ALL,
                    lonely_pairs=LonelyPairs.HEURISTIC,
                    energy_model="t04",
                    backend="baseopt",
                ),
            ),
            (linearfold, EnergyCfg(ctd=CtdCfg.D2, lonely_pairs=LonelyPairs.HEURISTIC)),
            (rnastructure, EnergyCfg(ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.HEURISTIC)),
            (viennarna, EnergyCfg(ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.HEURISTIC)),
            (viennarna, EnergyCfg(ctd=CtdCfg.D2, lonely_pairs=LonelyPairs.HEURISTIC)),
            (viennarna, EnergyCfg(ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.ON)),
            (viennarna, EnergyCfg(ctd=CtdCfg.D2, lonely_pairs=LonelyPairs.ON)),
            (sparsemfefold, EnergyCfg(ctd=CtdCfg.NONE, lonely_pairs=LonelyPairs.HEURISTIC)),
            (sparsernafold, EnergyCfg(ctd=CtdCfg.D2, lonely_pairs=LonelyPairs.HEURISTIC)),
        ]
        for program, _ in self.programs:
            program.limits.mem_bytes = mem_bytes_limit
            program.limits.time_sec = time_sec_limit

    def run(self) -> None:
        for program, cfg in self.programs:
            desc = program.desc(energy_cfg=cfg, subopt_cfg=None)
            dataset = self.memevault.dataset
            click.echo(f"Benchmarking folding with {desc} on {dataset}")
            output_path = self.output_dir / f"{dataset}_{desc}.results"
            if output_path.exists():
                raise RuntimeError(f"Output path {output_path} already exists")

            for rna_idx, rna in enumerate(self.memevault):
                click.echo(f"Running {program} on {rna_idx} {rna.name}")
                failed = False
                for run_idx in range(self.num_tries):
                    try:
                        _, res = program.fold(rna, cfg)
                    except Exception as e:
                        click.echo(f"Error running {program} on {rna.name}: {e}")
                        failed = True
                        break
                    df = pd.DataFrame(
                        {
                            "name": rna.name,
                            "run_idx": run_idx,
                            "length": len(rna),
                            "maxrss_bytes": res.maxrss_bytes,
                            "user_sec": res.user_sec,
                            "sys_sec": res.sys_sec,
                            "real_sec": res.real_sec,
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
