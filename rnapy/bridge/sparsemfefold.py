# Copyright 2022 Eliot Courtney.
from dataclasses import dataclass
from decimal import Decimal
from typing import override

from rnapy.bridge.rnapackage import RnaPackage
from rnapy.model.model_cfg import CtdCfg, EnergyCfg, LonelyPairs, SuboptCfg
from rnapy.model.parse.rna_parser import RnaParser
from rnapy.model.rna import Rna
from rnapy.util.command import CmdResult


@dataclass
class SparseMFEFold(RnaPackage):
    def _check_energy_cfg(self, cfg: EnergyCfg) -> None:
        if cfg.lonely_pairs != LonelyPairs.HEURISTIC:  # TODO(3): Check this.
            raise NotImplementedError(
                "SparseMFEFold does not support modifying lonely pairs behavior"
            )
        if cfg.ctd != CtdCfg.NONE:
            raise NotImplementedError("SparseMFEFold does not support turning on any CTDs")
        if cfg.energy_model is not None:
            raise NotImplementedError("SparseMFEFold energy model configuration not supported")

    @override
    def package_name(self) -> str:
        return "SparseMFEFold"

    @override
    def efn(self, rna: Rna, cfg: EnergyCfg) -> tuple[Decimal, CmdResult]:
        raise NotImplementedError

    @override
    def fold(self, rna: Rna, cfg: EnergyCfg) -> tuple[Rna, CmdResult]:
        self._check_energy_cfg(cfg)
        if rna.r is None:
            raise ValueError(f"RNA {rna.name} has no sequence")
        res = self._run_cmd("./src/SparseMFEFold", stdin_inp=rna.r, stdout_to_str=True)
        lines = res.stdout.strip().splitlines()
        if len(lines) < 2:
            raise ValueError(f"Unexpected SparseMFEFold output: {res.stdout!r}")
        seq, db_line = lines[0].strip(), lines[1].strip()
        db = db_line.split(maxsplit=1)[0]
        predicted = RnaParser.parse(name=rna.name, seq=seq, db=db)
        return predicted, res

    @override
    def partition(self, rna: Rna, cfg: EnergyCfg) -> None:
        raise NotImplementedError("SparseMFEFold does not support partition")

    @override
    def subopt(
        self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg
    ) -> tuple[list[Rna], CmdResult]:
        raise NotImplementedError("SparseMFEFold does not support suboptimal folding")
