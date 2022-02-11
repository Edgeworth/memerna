from dataclasses import dataclass

from scripts.bridge.rnapackage import RnaPackage
from scripts.model.config import EnergyCfg
from scripts.model.config import SuboptCfg
from scripts.model.parse import db_to_secondary
from scripts.model.rna import Rna


@dataclass
class MemeRna(RnaPackage):
    def energy_cfg_args(self, cfg: EnergyCfg):
        args = []
        if cfg.lonely_pairs:
            args.append("--lonely-pairs")
        else:
            args.append("--no-lonely-pairs")
        args += ["--ctd", f"{cfg.ctd}"]
        return args

    def subopt_cfg_args(self, cfg: SuboptCfg):
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

    def efn(self, rna: Rna, cfg: EnergyCfg):
        args = self.energy_cfg_args(cfg)
        res = self._run_cmd("./efn", *args, rna.r, rna.db())
        energy = float(res.stdout.splitlines()[0].strip()) / 10.0
        return energy, res

    def fold(self, rna: Rna, cfg: EnergyCfg):
        args = self.energy_cfg_args(cfg)
        res = self._run_cmd("./fold", *args, rna.r)
        _, db, _ = res.stdout.strip().split("\n")
        return Rna(rna.name, rna.r, db_to_secondary(db)), res

    def partition(self, rna: Rna, cfg: EnergyCfg):
        args = self.energy_cfg_args(cfg)
        raise NotImplementedError

    def subopt(self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg):
        args = self.energy_cfg_args(energy_cfg)
        args += self.subopt_cfg_args(subopt_cfg)
        res = self._run_cmd("./subopt", *args, rna.r)
        subopts = []
        for line in res.stdout.splitlines()[:-1]:
            energy, db = line.strip().split()
            subopts.append(Rna(name=rna.name, r=rna.r, s=db_to_secondary(db), energy=energy))
        return subopts, res

    def __str__(self):
        return "MemeRNA"
