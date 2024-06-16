# Copyright 2024 Eliot Courtney.
from dataclasses import dataclass
from decimal import Decimal
from typing import override

from rnapy.bridge.rnapackage import RnaPackage
from rnapy.model.model_cfg import EnergyCfg, LonelyPairs, SuboptCfg
from rnapy.model.parse.sequence import db_to_secondary
from rnapy.model.rna import Rna
from rnapy.util.command import CmdResult


@dataclass
class MemeRna01(RnaPackage):
    def _energy_cfg_args(self, cfg: EnergyCfg) -> list[str]:
        if cfg.lonely_pairs != LonelyPairs.HEURISTIC:
            raise NotImplementedError("memerna0.1 does not support modifying lonely pairs behavior")
        if cfg.energy_model is not None:
            raise NotImplementedError("memerna0.1 energy model configuration not supported")
        if cfg.ctd != "all":
            raise NotImplementedError("memerna0.1 only supports all CTDs")

        return []

    @override
    def package_name(self) -> str:
        return "memerna0.1"

    @override
    def efn(self, rna: Rna, cfg: EnergyCfg) -> tuple[Decimal, CmdResult]:
        raise NotImplementedError

    @override
    def fold(self, rna: Rna, cfg: EnergyCfg) -> tuple[Rna, CmdResult]:
        args = self._energy_cfg_args(cfg)
        if rna.r is None:
            raise ValueError(f"RNA {rna.name} has no sequence")
        res = self._run_cmd(
            "./build/fold",
            *args,
            "-memerna-data",
            str(self.path.joinpath("data")),
            rna.r,
            stdout_to_str=True,
        )
        _, db, _ = res.stdout.strip().split("\n")
        return Rna(rna.name, rna.r, db_to_secondary(db)), res

    @override
    def partition(self, rna: Rna, cfg: EnergyCfg) -> None:
        raise NotImplementedError

    @override
    def subopt(
        self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg
    ) -> tuple[list[Rna], CmdResult]:
        raise NotImplementedError
