from pathlib import Path
from rnapy.bridge.rnapackage import RnaPackage
from rnapy.data.memevault import MemeVault
from rnapy.bridge.memerna import MemeRna
from rnapy.bridge.rnastructure import RNAstructure
from rnapy.bridge.sparsemfefold import SparseMfeFold
from rnapy.bridge.viennarna import ViennaRna
import click

from rnapy.model.model_cfg import CtdCfg, EnergyCfg


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
        rnastructure: RNAstructure,
        viennarna: ViennaRna,
        sparsemfefold: SparseMfeFold,
    ) -> None:
        self.memevault = memevault
        self.output_dir = output_dir
        self.programs = [
            (memerna, EnergyCfg(), memerna.name()),
            (rnastructure, EnergyCfg(), rnastructure.name()),
            (viennarna, EnergyCfg(), viennarna.name() + "-d3"),
            (viennarna, EnergyCfg(ctd=CtdCfg.D2), viennarna.name() + "-d2"),
            (sparsemfefold, EnergyCfg(), sparsemfefold.name()),
        ]

    def run(self) -> None:
        for program, cfg, name in self.programs:
            dataset = self.memevault.dataset
            click.echo(f"Benchmarking folding with {name} on {dataset}")
            output_path = self.output_dir / f"{name}_{dataset}_fold.results"
            assert not output_path.exists(), f"Output path {output_path} already exists"

            with open(output_path, "w", encoding="utf-8") as f:
                for rna_idx, rna in enumerate(self.memevault):
                    click.echo(f"Running {program} on {rna_idx} {rna.name}")
                    for run_idx in range(self.BENCHMARK_NUM_TRIES):
                        _, res = program.fold(rna, cfg)
                        f.write(
                            f"{rna.name} {run_idx} {len(rna)} {res.maxrss_bytes} "
                            f"{res.user_sec} {res.sys_sec} {res.real_sec}\n"
                        )
