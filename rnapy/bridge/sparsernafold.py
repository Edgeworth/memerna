# Copyright 2023 Eliot Courtney.
from dataclasses import dataclass
from decimal import Decimal
from typing import override

from rnapy.bridge.rnapackage import RnaPackage
from rnapy.model.model_cfg import CtdCfg, EnergyCfg, LonelyPairs, SuboptCfg
from rnapy.model.parse.rna_parser import RnaParser
from rnapy.model.rna import Rna
from rnapy.util.command import CmdResult


@dataclass
class SparseRNAFolD(RnaPackage):
    def _energy_cfg_args(self, cfg: EnergyCfg) -> list[str]:
        args = []

        if cfg.lonely_pairs != LonelyPairs.HEURISTIC:
            raise NotImplementedError(
                "SparseRNAFolD does not support modifying lonely pairs behavior"
            )

        match cfg.ctd:
            case CtdCfg.NONE:
                args.append("-d0")
            case CtdCfg.D2:
                args.append("-d2")
            case CtdCfg.NO_COAX:
                raise NotImplementedError(
                    "SparseRNAFolD does not support CTDs with no coaxial stacking"
                )
            case CtdCfg.ALL:
                raise NotImplementedError("SparseRNAFolD does not support all CTDs")

        if cfg.model is not None:
            raise NotImplementedError("SparseRNAFolD energy model configuration not supported")
        return args

    @override
    def name(self) -> str:
        return "SparseRNAFolD"

    @override
    def efn(self, rna: Rna, cfg: EnergyCfg) -> tuple[Decimal, CmdResult]:
        raise NotImplementedError

    @override
    def fold(self, rna: Rna, cfg: EnergyCfg) -> tuple[Rna, CmdResult]:
        args = self._energy_cfg_args(cfg)
        if rna.r is None:
            raise ValueError(f"RNA {rna.name} has no sequence")
        res = self._run_cmd("./build/SparseRNAFolD", *args, rna.r, stdout_to_str=True)
        seq, db = res.stdout.strip().split("\n")
        db = db.split(" ")[0]
        predicted = RnaParser.parse(name=rna.name, seq=seq.strip(), db=db.strip())
        return predicted, res

    @override
    def partition(self, rna: Rna, cfg: EnergyCfg) -> None:
        raise NotImplementedError

    @override
    def subopt(
        self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg
    ) -> tuple[list[Rna], CmdResult]:
        raise NotImplementedError

    def __str__(self) -> str:
        return self.name()
