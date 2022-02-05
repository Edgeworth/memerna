# Copyright 2022 Eliot Courtney.

from scripts.model.config import EnergyCfg, SuboptCfg
from scripts.model.rna import Rna


class RnaPackage:
    def efn(self, rna: Rna):
        pass

    def fold(self, rna: Rna, cfg: EnergyCfg):
        pass

    def subopt(self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg):
        pass
