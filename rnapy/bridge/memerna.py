# Copyright 2022 Eliot Courtney.
from dataclasses import dataclass
from decimal import Decimal

from rnapy.bridge.rnapackage import RnaPackage
from rnapy.model.model_cfg import EnergyCfg, SuboptCfg
from rnapy.model.parse.sequence import db_to_secondary
from rnapy.model.rna import Rna
from rnapy.util.command import CmdResult


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
        args += ["--no-ctd-output"]
        return args

    def name(self) -> str:
        return "memerna"

    def efn(self, rna: Rna, cfg: EnergyCfg) -> tuple[Decimal, CmdResult]:
        args = self._energy_cfg_args(cfg)
        if rna.r is None:
            raise ValueError(f"RNA {rna.name} has no sequence")
        res = self._run_cmd("./efn", *args, rna.r, rna.db())
        energy = Decimal(res.stdout.splitlines()[0].strip())
        return energy, res

    def fold(self, rna: Rna, cfg: EnergyCfg) -> tuple[Rna, CmdResult]:
        args = self._energy_cfg_args(cfg)
        if rna.r is None:
            raise ValueError(f"RNA {rna.name} has no sequence")
        res = self._run_cmd("./fold", *args, rna.r)
        _, db, _ = res.stdout.strip().split("\n")
        return Rna(rna.name, rna.r, db_to_secondary(db)), res

    def partition(self, rna: Rna, cfg: EnergyCfg) -> None:
        raise NotImplementedError

    def subopt(
        self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg
    ) -> tuple[list[Rna], CmdResult]:
        args = self._energy_cfg_args(energy_cfg)
        args += self._subopt_cfg_args(subopt_cfg)
        if rna.r is None:
            raise ValueError(f"RNA {rna.name} has no sequence")
        res = self._run_cmd("./subopt", *args, rna.r)
        subopts = []
        for line in res.stdout.splitlines()[:-1]:
            energy_str, db = line.strip().split()
            energy = Decimal(energy_str)
            subopts.append(Rna(name=rna.name, r=rna.r, s=db_to_secondary(db), energy=energy))
        return subopts, res

    def __str__(self) -> str:
        return "memerna"
