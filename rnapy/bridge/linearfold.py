# Copyright 2023 E.
from dataclasses import dataclass
from decimal import Decimal

from rnapy.bridge.rnapackage import RnaPackage
from rnapy.model.model_cfg import CtdCfg
from rnapy.model.model_cfg import EnergyCfg
from rnapy.model.model_cfg import LonelyPairs
from rnapy.model.model_cfg import SuboptCfg
from rnapy.model.parse.rna_parser import RnaParser
from rnapy.model.rna import Rna
from rnapy.util.command import CmdResult


@dataclass
class LinearFold(RnaPackage):
    def _energy_cfg_args(self, cfg: EnergyCfg) -> list[str]:
        args = []

        if cfg.lonely_pairs != LonelyPairs.HEURISTIC:
            raise NotImplementedError("LinearFold does not support turning on lonely pairs")

        match cfg.ctd:
            case CtdCfg.NONE:
                args.append("-d0")
            case CtdCfg.D2:
                args.append("-d2")
            case CtdCfg.NO_COAX:
                raise NotImplementedError(
                    "LinearFold does not support CTDs with no coaxial stacking",
                )
            case CtdCfg.ALL:
                raise NotImplementedError("LinearFold does not support all CTDs")

        if cfg.model is not None:
            raise NotImplementedError("LinearFold energy model configuration not supported")
        return args

    def name(self) -> str:
        return "LinearFold"

    def efn(self, rna: Rna, cfg: EnergyCfg) -> tuple[Decimal, CmdResult]:
        raise NotImplementedError

    def fold(self, rna: Rna, cfg: EnergyCfg) -> tuple[Rna, CmdResult]:
        args = self._energy_cfg_args(cfg)
        res = self._run_cmd("./linearfold", *args, "-V", inp=rna.r)
        seq, db = res.stdout.strip().split("\n")
        db = db.split(" ")[0]
        predicted = RnaParser.parse(name=rna.name, seq=seq.strip(), db=db.strip())
        return predicted, res

    def partition(self, rna: Rna, cfg: EnergyCfg) -> None:
        raise NotImplementedError

    def subopt(
        self,
        rna: Rna,
        energy_cfg: EnergyCfg,
        subopt_cfg: SuboptCfg,
    ) -> tuple[list[Rna], CmdResult]:
        raise NotImplementedError

    def __str__(self) -> str:
        return "LinearFold"
