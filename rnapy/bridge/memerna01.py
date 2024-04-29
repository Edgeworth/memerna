# Copyright 2024 Eliot Courtney.
from dataclasses import dataclass
from decimal import Decimal

from rnapy.bridge.rnapackage import RnaPackage
from rnapy.model.model_cfg import EnergyCfg, LonelyPairs, SuboptCfg
from rnapy.model.parse.sequence import db_to_secondary
from rnapy.model.rna import Rna
from rnapy.util.command import CmdResult


@dataclass
class MemeRna01(RnaPackage):
    def _energy_cfg_args(self, cfg: EnergyCfg) -> list[str]:
        if cfg.lonely_pairs != LonelyPairs.HEURISTIC:
            raise NotImplementedError("memerna0.1 does not support turning on lonely pairs")
        if cfg.model is not None:
            raise NotImplementedError("memerna0.1 energy model configuration not supported")
        if cfg.ctd != "all":
            raise NotImplementedError("memerna0.1 only supports all CTDs")

        return []

    def name(self) -> str:
        return "memerna"

    def efn(self, rna: Rna, cfg: EnergyCfg) -> tuple[Decimal, CmdResult]:
        raise NotImplementedError

    def fold(self, rna: Rna, cfg: EnergyCfg) -> tuple[Rna, CmdResult]:
        args = self._energy_cfg_args(cfg)
        if rna.r is None:
            raise ValueError(f"RNA {rna.name} has no sequence")
        res = self._run_cmd(
            "./build/fold", *args, "-memerna-data", str(self.path.joinpath("data")), rna.r
        )
        _, db, _ = res.stdout.strip().split("\n")
        return Rna(rna.name, rna.r, db_to_secondary(db)), res

    def partition(self, rna: Rna, cfg: EnergyCfg) -> None:
        raise NotImplementedError

    def subopt(
        self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg
    ) -> tuple[list[Rna], CmdResult]:
        raise NotImplementedError

    def __str__(self) -> str:
        return "memerna"
