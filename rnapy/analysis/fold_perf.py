# Copyright 2022 Eliot Courtney.
from pathlib import Path
from typing import Any

import click
import polars as pl

from rnapy.bridge.linearfold import LinearFold
from rnapy.bridge.memerna01 import MemeRna01
from rnapy.bridge.rnapackage import RnaPackage
from rnapy.bridge.rnastructure import RNAstructure
from rnapy.bridge.sparsemfefold import SparseMFEFold
from rnapy.bridge.sparsernafold import SparseRNAFolD
from rnapy.bridge.viennarna import ViennaRna
from rnapy.data.memevault import MemeVault
from rnapy.model.model_cfg import CtdCfg, EnergyCfg, LonelyPairs
from rnapy.model.rna import Rna
from rnapy.util.util import append_ndjson, row_by_key, strict_merge


class FoldPerfRunner:
    num_tries: int
    memevault: MemeVault
    output_path: Path
    programs: list[tuple[RnaPackage, EnergyCfg]]

    def __init__(
        self,
        *,
        num_tries: int,
        memevault: MemeVault,
        output_path: Path,
        memerna01: MemeRna01,
        linearfold: LinearFold,
        rnastructure: RNAstructure,
        viennarna: ViennaRna,
        sparsemfefold: SparseMFEFold,
        sparsernafold: SparseRNAFolD,
    ) -> None:
        self.num_tries = num_tries
        self.memevault = memevault
        self.output_path = output_path
        self.programs = [
            (memerna01, EnergyCfg(ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.HEURISTIC)),
            (linearfold, EnergyCfg(ctd=CtdCfg.D2, lonely_pairs=LonelyPairs.HEURISTIC)),
            (rnastructure, EnergyCfg(ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.HEURISTIC)),
            (viennarna, EnergyCfg(ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.HEURISTIC)),
            (viennarna, EnergyCfg(ctd=CtdCfg.D2, lonely_pairs=LonelyPairs.HEURISTIC)),
            (viennarna, EnergyCfg(ctd=CtdCfg.ALL, lonely_pairs=LonelyPairs.ON)),
            (viennarna, EnergyCfg(ctd=CtdCfg.D2, lonely_pairs=LonelyPairs.ON)),
            (sparsemfefold, EnergyCfg(ctd=CtdCfg.NONE, lonely_pairs=LonelyPairs.HEURISTIC)),
            (sparsernafold, EnergyCfg(ctd=CtdCfg.D2, lonely_pairs=LonelyPairs.HEURISTIC)),
        ]

    def _run_once(self, program: RnaPackage, cfg: EnergyCfg, rna: Rna) -> None:
        dataset = self.memevault.dataset
        desc = program.desc(energy_cfg=cfg, subopt_cfg=None)
        name = rna.require_name()

        for run_idx in range(self.num_tries):
            data_keys = strict_merge(
                desc,
                {"dataset": dataset, "rna_name": name, "rna_length": len(rna), "run_idx": run_idx},
            )

            row = row_by_key(self.output_path, data_keys)
            if row is not None:
                if row["failed"]:
                    click.echo(f"Skipping run {row} as it failed previously.")
                    return
                click.echo(f"Skipping run {row} as it already exists in {self.output_path}")
                continue

            failed = False
            data_values: dict[str, Any] = {}
            try:
                folded, res = program.fold(rna, cfg)
                data_values = {
                    "energy": folded.energy,
                    "maxrss_bytes": res.maxrss_bytes,
                    "user_sec": res.user_sec,
                    "sys_sec": res.sys_sec,
                    "real_sec": res.real_sec,
                }
            except NotImplementedError:
                # Packages raise this for unsupported configurations, which the caller skips.
                raise
            except Exception as e:
                click.echo(f"Error running {program} on {name}: {e}")
                failed = True

            data = strict_merge(data_keys, data_values, {"failed": failed})
            append_ndjson(self.output_path, pl.DataFrame([data]))

            if failed:
                return

    def run(self) -> None:
        for program, cfg in self.programs:
            desc = program.desc(energy_cfg=cfg, subopt_cfg=None)
            dataset = self.memevault.dataset
            click.echo(f"Benchmarking folding with {desc} on {dataset}")
            try:
                for rna_idx, rna in enumerate(self.memevault):
                    click.echo(f"Running {program} on {rna_idx} {rna.name}")
                    self._run_once(program, cfg, rna)
            except NotImplementedError as e:
                click.echo(f"Skipping {desc}: {e}")
