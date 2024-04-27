# Copyright 2022 Eliot Courtney.
from pathlib import Path

import click
import pandas as pd

from rnapy.bridge.linearfold import LinearFold
from rnapy.bridge.memerna import MemeRna
from rnapy.bridge.rnapackage import RnaPackage
from rnapy.bridge.rnastructure import RNAstructure
from rnapy.bridge.sparsemfefold import SparseMfeFold
from rnapy.bridge.sparsernafold import SparseRNAFolD
from rnapy.bridge.viennarna import ViennaRna
from rnapy.data.memevault import MemeVault
from rnapy.model.model_cfg import CtdCfg, EnergyCfg, LonelyPairs


class FoldPerfRunner:
    num_tries: int
    memevault: MemeVault
    output_dir: Path
    programs: list[tuple[RnaPackage, EnergyCfg, str]]

    def __init__(
        self,
        *,
        num_tries: int,
        time_sec_limit: int | None,
        mem_bytes_limit: int | None,
        memevault: MemeVault,
        output_dir: Path,
        memerna: MemeRna,
        linearfold: LinearFold,
        rnastructure: RNAstructure,
        viennarna: ViennaRna,
        sparsemfefold: SparseMfeFold,
        sparsernafold: SparseRNAFolD,
    ) -> None:
        self.num_tries = num_tries
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
            (sparsernafold, EnergyCfg(ctd=CtdCfg.D2), sparsernafold.name()),
        ]
        for program, _, _ in self.programs:
            program.limits.mem_bytes = mem_bytes_limit
            program.limits.time_sec = time_sec_limit

    def run(self) -> None:
        for program, cfg, name in self.programs:
            dataset = self.memevault.dataset
            click.echo(f"Benchmarking folding with {name} on {dataset}")
            output_path = self.output_dir / f"{dataset}_{name}.results"
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
