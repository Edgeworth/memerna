# Copyright 2022 Eliot Courtney.
from dataclasses import dataclass

from scripts.bridge.rnapackage import RnaPackage
from scripts.model.config import CtdCfg
from scripts.model.config import EnergyCfg
from scripts.model.config import SuboptCfg
from scripts.model.rna import Rna
from scripts.util.command import CmdResult


@dataclass
class SparseMfeFold(RnaPackage):
    def check_energy_cfg(self, cfg: EnergyCfg) -> None:
        if cfg.lonely_pairs:  # TODO: Check this.
            raise NotImplementedError("SparseMFEFold does not support turning on lonely pairs")
        if cfg.ctd != CtdCfg.NONE:
            raise NotImplementedError("SparseMFEFold does not support turning on any CTDs")

    def efn(self, rna: Rna, cfg: EnergyCfg) -> tuple[float, CmdResult]:
        raise NotImplementedError

    def fold(self, rna: Rna, cfg: EnergyCfg) -> tuple[Rna, CmdResult]:
        self.check_energy_cfg(cfg)
        res = self._run_cmd("./src/SparseMFEFold", inp=rna.r)
        seq, db = res.stdout.strip().split("\n")
        db = db.split(" ")[0]
        predicted = Rna.parse(name=rna.name, seq=seq.strip(), db=db.strip())
        return predicted, res

    def partition(self, rna: Rna, cfg: EnergyCfg) -> None:
        raise NotImplementedError("SparseMFEFold does not support partition")

    def subopt(
        self,
        rna: Rna,
        energy_cfg: EnergyCfg,
        subopt_cfg: SuboptCfg,
    ) -> tuple[list[Rna], CmdResult]:
        raise NotImplementedError("SparseMFEFold does not support suboptimal folding")

    def __str__(self) -> str:
        return "SparseMFEFold"
