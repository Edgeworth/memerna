# Copyright 2023 Eliot Courtney.
from pathlib import Path

import click
import pandas as pd
from rnapy.analysis.metrics import RnaAccuracy
from rnapy.bridge.memerna import MemeRna
from rnapy.bridge.rnapackage import RnaPackage
from rnapy.bridge.rnastructure import RNAstructure
from rnapy.bridge.sparsemfefold import SparseMfeFold
from rnapy.bridge.viennarna import ViennaRna
from rnapy.data.memevault import MemeVault
from rnapy.model.model_cfg import CtdCfg
from rnapy.model.model_cfg import EnergyCfg


class FoldAccuracyRunner:
    memevault: MemeVault
    output_dir: Path
    programs: list[tuple[RnaPackage, EnergyCfg, str]]

    def __init__(
        self,
        memevault: MemeVault,
        output_dir: Path,
        memerna: MemeRna,
    ) -> None:
        self.memevault = memevault
        self.output_dir = output_dir
        self.programs = [
            # (memerna, EnergyCfg(model="t04p1"), memerna.name() + "-t04p1"),
            (
                memerna,
                EnergyCfg(model="t04p2", lonely_pairs=True),
                memerna.name() + "-t04p2-lonely",
            ),
            # (memerna, EnergyCfg(model="t04p2"), memerna.name() + "-t04p2"),
            # (memerna, EnergyCfg(model="t22p2"), memerna.name() + "-t22p2"),
        ]

    @staticmethod
    def extract_family(name: str) -> str:
        return name.split("_")[0]

    def run(self) -> None:
        for program, cfg, name in self.programs:
            dataset = self.memevault.dataset
            click.echo(f"Folding accuracy with {name} on {dataset}")
            output_path = self.output_dir / f"{dataset}_{name}.results"
            assert not output_path.exists(), f"Output path {output_path} already exists"

            for rna_idx, rna in enumerate(self.memevault):
                click.echo(f"Running {name} on {rna_idx} {rna.name}")
                assert rna.name
                pred, _ = program.fold(rna, cfg)
                accuracy = RnaAccuracy.from_rna(rna, pred)
                df = pd.DataFrame(
                    {
                        "name": rna.name,
                        "family": self.extract_family(rna.name),
                        "ppv": accuracy.ppv,
                        "sensitivity": accuracy.sensitivity,
                        "f1": accuracy.f1_score(),
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
