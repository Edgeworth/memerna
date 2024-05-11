# Copyright 2022 Eliot Courtney.
import tempfile
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import override

from rnapy.bridge.rnapackage import RnaPackage
from rnapy.model.model_cfg import EnergyCfg, SuboptCfg
from rnapy.model.parse.sequence import db_to_secondary
from rnapy.model.rna import Rna
from rnapy.util.command import CmdResult
from rnapy.util.util import fast_linecount


@dataclass
class MemeRna(RnaPackage):
    def _energy_cfg_args(self, cfg: EnergyCfg) -> list[str]:
        args = []

        if cfg.model:
            args += ["--energy-model", f"{cfg.model}"]
        args += ["--ctd", f"{cfg.ctd}"]
        args += ["--lonely-pairs", f"{cfg.lonely_pairs}"]
        return args

    def _subopt_cfg_args(self, cfg: SuboptCfg) -> list[str]:
        args = []
        if cfg.delta:
            args += ["--subopt-delta", str(cfg.delta)]
        if cfg.strucs:
            args += ["--subopt-strucs", str(cfg.strucs)]
        if cfg.time_secs:
            args += ["--subopt-time-secs", str(cfg.time_secs)]
        if cfg.sorted_strucs:
            args += ["--subopt-sorted"]
        else:
            args += ["--no-subopt-sorted"]
        if cfg.algorithm:
            args += ["--subopt-alg", f"{cfg.algorithm}"]
        args += ["--no-ctd-output"]
        return args

    @override
    def name(self) -> str:
        return "memerna"

    @override
    def efn(self, rna: Rna, cfg: EnergyCfg) -> tuple[Decimal, CmdResult]:
        args = self._energy_cfg_args(cfg)
        if rna.r is None:
            raise ValueError(f"RNA {rna.name} has no sequence")
        res = self._run_cmd("./efn", *args, rna.r, rna.db(), stdout_to_str=True)
        energy = Decimal(res.stdout.splitlines()[0].strip())
        return energy, res

    @override
    def fold(self, rna: Rna, cfg: EnergyCfg) -> tuple[Rna, CmdResult]:
        args = self._energy_cfg_args(cfg)
        if rna.r is None:
            raise ValueError(f"RNA {rna.name} has no sequence")
        res = self._run_cmd("./fold", *args, rna.r, stdout_to_str=True)
        _, db, _ = res.stdout.strip().split("\n")
        return Rna(rna.name, rna.r, db_to_secondary(db)), res

    @override
    def partition(self, rna: Rna, cfg: EnergyCfg) -> None:
        raise NotImplementedError

    @override
    def subopt(
        self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg
    ) -> tuple[list[Rna] | int, CmdResult]:
        args = self._energy_cfg_args(energy_cfg)
        args += self._subopt_cfg_args(subopt_cfg)
        if rna.r is None:
            raise ValueError(f"RNA {rna.name} has no sequence")

        with tempfile.NamedTemporaryFile("w") as f:
            stdout_path = Path(f.name)
            res = self._run_cmd(
                "./subopt", *args, rna.r, stdout_to_str=False, stdout_path=stdout_path
            )
            if subopt_cfg.count_only:
                count = fast_linecount(stdout_path) - 1
                return count, res

            subopts = []
            for line in stdout_path.read_text().splitlines()[:-1]:
                energy_str, db = line.strip().split()
                energy = Decimal(energy_str)
                subopts.append(Rna(name=rna.name, r=rna.r, s=db_to_secondary(db), energy=energy))
            return subopts, res

    def __str__(self) -> str:
        return self.name()
