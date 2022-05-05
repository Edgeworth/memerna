# Copyright 2022 Eliot Courtney.
from dataclasses import dataclass

from scripts.bridge.rnapackage import RnaPackage
from scripts.model.config import EnergyCfg
from scripts.model.config import SuboptCfg
from scripts.model.parse.sequence import db_to_secondary
from scripts.model.rna import Rna
from scripts.util.command import CmdResult


@dataclass
class MemeRna(RnaPackage):
    def energy_cfg_args(self, cfg: EnergyCfg) -> list[str]:
        args = []
        if cfg.lonely_pairs:
            args.append("--lonely-pairs")
        else:
            args.append("--no-lonely-pairs")
        args += ["--ctd", f"{cfg.ctd}"]
        return args

    def subopt_cfg_args(self, cfg: SuboptCfg) -> list[str]:
        args = []
        if cfg.delta:
            args += ["--subopt-delta", str(cfg.delta)]
        if cfg.strucs:
            args += ["--subopt-strucs", str(cfg.strucs)]
        if cfg.sorted:
            args += ["--subopt-sorted"]
        else:
            args += ["--no-subopt-sorted"]
        return args

    def efn(self, rna: Rna, cfg: EnergyCfg) -> tuple[float, CmdResult]:
        args = self.energy_cfg_args(cfg)
        assert rna.r is not None
        res = self._run_cmd("./efn", *args, rna.r, rna.db())
        energy = float(res.stdout.splitlines()[0].strip()) / 10.0
        return energy, res

    def fold(self, rna: Rna, cfg: EnergyCfg) -> tuple[Rna, CmdResult]:
        args = self.energy_cfg_args(cfg)
        assert rna.r is not None
        res = self._run_cmd("./fold", *args, rna.r)
        _, db, _ = res.stdout.strip().split("\n")
        return Rna(rna.name, rna.r, db_to_secondary(db)), res

    def partition(self, rna: Rna, cfg: EnergyCfg) -> None:
        raise NotImplementedError

    def subopt(
        self,
        rna: Rna,
        energy_cfg: EnergyCfg,
        subopt_cfg: SuboptCfg,
    ) -> tuple[list[Rna], CmdResult]:
        args = self.energy_cfg_args(energy_cfg)
        args += self.subopt_cfg_args(subopt_cfg)
        assert rna.r is not None
        res = self._run_cmd("./subopt", *args, rna.r)
        subopts = []
        for line in res.stdout.splitlines()[:-1]:
            energy_str, db = line.strip().split()
            energy = int(energy_str)
            subopts.append(Rna(name=rna.name, r=rna.r, s=db_to_secondary(db), energy=energy))
        return subopts, res

    def __str__(self) -> str:
        return "MemeRNA"
