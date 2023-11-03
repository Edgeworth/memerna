# Copyright 2022 Eliot Courtney.
from pathlib import Path

import click
import pandas as pd

from rnapy.bridge.linearfold import LinearFold
from rnapy.bridge.memerna import MemeRna
from rnapy.bridge.rnapackage import RnaPackage
from rnapy.bridge.rnastructure import RNAstructure
from rnapy.bridge.sparsemfefold import SparseMfeFold
from rnapy.bridge.viennarna import ViennaRna
from rnapy.data.memevault import MemeVault
from rnapy.model.model_cfg import CtdCfg, EnergyCfg, LonelyPairs


class FoldPerfRunner:
    BENCHMARK_NUM_TRIES = 5

    memevault: MemeVault
    output_dir: Path
    programs: list[tuple[RnaPackage, EnergyCfg, str]]

    def __init__(
        self,
        memevault: MemeVault,
        output_dir: Path,
        memerna: MemeRna,
        linearfold: LinearFold,
        rnastructure: RNAstructure,
        viennarna: ViennaRna,
        sparsemfefold: SparseMfeFold,
    ) -> None:
        self.memevault = memevault
        self.output_dir = output_dir
        self.programs = [
            (memerna, EnergyCfg(), memerna.name()),
            (linearfold, EnergyCfg(ctd=CtdCfg.D2), linearfold.name()),
            (rnastructure, EnergyCfg(), rnastructure.name()),
            (viennarna, EnergyCfg(), viennarna.name() + "-d3-noLP"),
            (viennarna, EnergyCfg(ctd=CtdCfg.D2), viennarna.name() + "-d2-noLP"),
            (viennarna, EnergyCfg(lonely_pairs=LonelyPairs.HEURISTIC), viennarna.name() + "-d3"),
            (
                viennarna,
                EnergyCfg(lonely_pairs=LonelyPairs.HEURISTIC, ctd=CtdCfg.D2),
                viennarna.name() + "-d2",
            ),
            (sparsemfefold, EnergyCfg(ctd=CtdCfg.NONE), sparsemfefold.name()),
        ]

    def run(self) -> None:
        for program, cfg, name in self.programs:
            dataset = self.memevault.dataset
            click.echo(f"Benchmarking folding with {name} on {dataset}")
            output_path = self.output_dir / f"{dataset}_{name}.results"
            assert not output_path.exists(), f"Output path {output_path} already exists"

            for rna_idx, rna in enumerate(self.memevault):
                click.echo(f"Running {program} on {rna_idx} {rna.name}")
                for run_idx in range(self.BENCHMARK_NUM_TRIES):
                    _, res = program.fold(rna, cfg)
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
