# Copyright 2023 Eliot Courtney.
from pathlib import Path

import click
import polars as pl

from rnapy.analysis.metrics import RnaAccuracy
from rnapy.bridge.memerna import MemeRna
from rnapy.bridge.rnapackage import RnaPackage
from rnapy.data.memevault import MemeVault
from rnapy.model.model_cfg import CtdCfg, EnergyCfg, LonelyPairs
from rnapy.util.util import append_csv


class FoldAccuracyRunner:
    memevault: MemeVault
    output_dir: Path
    programs: list[tuple[RnaPackage, EnergyCfg]]

    def __init__(self, memevault: MemeVault, output_dir: Path, memerna: MemeRna) -> None:
        self.memevault = memevault
        self.output_dir = output_dir
        self.programs = [
            (
                memerna,
                EnergyCfg(
                    ctd=CtdCfg.ALL,
                    lonely_pairs=LonelyPairs.HEURISTIC,
                    energy_model="t04",
                    backend="base",
                ),
            ),
            (
                memerna,
                EnergyCfg(
                    ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.ON, energy_model="t04", backend="base"
                ),
            ),
            (
                memerna,
                EnergyCfg(
                    ctd=CtdCfg.NO_COAX,
                    lonely_pairs=LonelyPairs.HEURISTIC,
                    energy_model="t04",
                    backend="base",
                ),
            ),
            (
                memerna,
                EnergyCfg(
                    ctd=CtdCfg.ALL,
                    lonely_pairs=LonelyPairs.HEURISTIC,
                    energy_model="t22",
                    backend="stack",
                ),
            ),
        ]

    @staticmethod
    def extract_family(name: str) -> str:
        return name.split("_")[0]

    def run(self) -> None:
        for program, cfg in self.programs:
            dataset = self.memevault.dataset
            desc = program.desc(energy_cfg=cfg, subopt_cfg=None)
            click.echo(f"Folding accuracy with {desc} on {dataset}")
            output_path = self.output_dir / f"{dataset}_{desc}.results"
            if output_path.exists():
                raise RuntimeError(f"Output path {output_path} already exists")

            for rna_idx, rna in enumerate(self.memevault):
                click.echo(f"Running {desc} on {rna_idx} {rna.name}")
                if not rna.name:
                    raise ValueError(f"RNA name is empty: {rna}")
                pred, _ = program.fold(rna, cfg)
                accuracy = RnaAccuracy.from_rna(rna, pred)
                df = pl.DataFrame(
                    {
                        "name": [rna.name],
                        "family": [self.extract_family(rna.name)],
                        "ppv": [accuracy.ppv],
                        "sensitivity": [accuracy.sensitivity],
                        "f1": [accuracy.f1_score()],
                    }
                )
                append_csv(output_path, df)
