# Copyright 2023 Eliot Courtney.
from pathlib import Path

import click
import polars as pl

from rnapy.analysis.metrics import RnaAccuracy
from rnapy.bridge.memerna import MemeRna
from rnapy.bridge.rnapackage import RnaPackage
from rnapy.data.memevault import MemeVault
from rnapy.model.model_cfg import CtdCfg, EnergyCfg, LonelyPairs
from rnapy.model.rna import Rna
from rnapy.util.util import append_ndjson, row_by_key, strict_merge


class FoldAccuracyRunner:
    memevault: MemeVault
    output_path: Path
    programs: list[tuple[RnaPackage, EnergyCfg]]

    def __init__(self, *, memevault: MemeVault, output_path: Path, memerna: MemeRna) -> None:
        self.memevault = memevault
        self.output_path = output_path
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
    def _extract_family(name: str) -> str:
        return name.split("_")[0]

    def _run_once(self, program: RnaPackage, cfg: EnergyCfg, rna: Rna) -> None:
        desc = program.desc(energy_cfg=cfg, subopt_cfg=None)
        dataset = self.memevault.dataset

        if not rna.name:
            raise ValueError(f"RNA name is empty: {rna}")

        data_keys = strict_merge(
            desc, {"dataset": dataset, "rna_name": rna.name, "rna_length": len(rna)}
        )

        row = row_by_key(self.output_path, data_keys)
        if row is not None:
            click.echo(f"Skipping run {row} as it already exists in {self.output_path}")
            return

        failed = False
        data_values: dict = {}
        try:
            pred, _ = program.fold(rna, cfg)
            accuracy = RnaAccuracy.from_rna(rna, pred)
            data_values = {
                "family": self._extract_family(rna.name),
                "ppv": accuracy.ppv,
                "sensitivity": accuracy.sensitivity,
                "f1": accuracy.f1_score(),
            }
        except Exception as e:
            click.echo(f"Error running {program} on {rna.name}: {e}")
            failed = True

        data = strict_merge(data_keys, data_values, {"failed": failed})
        append_ndjson(self.output_path, pl.DataFrame([data]))

    def run(self) -> None:
        for program, cfg in self.programs:
            dataset = self.memevault.dataset
            desc = program.desc(energy_cfg=cfg, subopt_cfg=None)
            click.echo(f"Folding accuracy with {desc} on {dataset}")
            for rna_idx, rna in enumerate(self.memevault):
                click.echo(f"Running {desc} on {rna_idx} {rna.name}")
                self._run_once(program, cfg, rna)
