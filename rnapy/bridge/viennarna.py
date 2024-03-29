# Copyright 2022 E.
import re
import tempfile
from dataclasses import dataclass
from decimal import Decimal

from rnapy.bridge.rnapackage import RnaPackage
from rnapy.model.model_cfg import CtdCfg, EnergyCfg, LonelyPairs, SuboptCfg
from rnapy.model.parse.rna_parser import RnaParser
from rnapy.model.parse.sequence import db_to_secondary
from rnapy.model.rna import Rna
from rnapy.util.command import CmdResult


@dataclass
class ViennaRna(RnaPackage):
    def _energy_cfg_args(self, cfg: EnergyCfg) -> list[str]:
        args = []

        match cfg.lonely_pairs:
            case LonelyPairs.OFF:
                args.append("--noLP")
            case LonelyPairs.HEURISTIC:
                pass
            case LonelyPairs.ON:
                raise NotImplementedError(
                    "ViennaRNA does not support non-heuristic disallowing of lonely pairs"
                )
        match cfg.ctd:
            case CtdCfg.NONE:
                args.append("-d0")
            case CtdCfg.D2:
                args.append("-d2")
            case CtdCfg.NO_COAX:
                raise NotImplementedError(
                    "ViennaRNA does not support CTDs with no coaxial stacking"
                )
            case CtdCfg.ALL:
                args.append("-d3")
        if cfg.model is not None:
            raise NotImplementedError("ViennaRNA energy model configuration not supported")
        return args

    def _subopt_cfg_args(self, cfg: SuboptCfg) -> list[str]:
        args = []
        if cfg.delta:
            args += ["--deltaEnergy", f"{cfg.delta}"]
        if cfg.strucs:
            raise NotImplementedError(
                "ViennaRNA does not support reporting a maximum number of suboptimal structures"
            )
        if cfg.time_secs:
            raise NotImplementedError("ViennaRNA does not support reporting a maximum running time")
        if cfg.sorted:
            args += ["--sorted"]
        return args

    def name(self) -> str:
        return "ViennaRNA"

    def efn(self, rna: Rna, cfg: EnergyCfg) -> tuple[Decimal, CmdResult]:
        args = self._energy_cfg_args(cfg)
        res = self._run_cmd("./src/bin/RNAeval", *args, inp=f"{rna.r}\n{rna.db()}")
        match = re.search(r"\s+\(\s*([0-9\.\-]+)\s*\)", res.stdout.strip())
        if match is None:
            raise ValueError(f"Could not find energy in {res.stdout}")
        energy = Decimal(match.group(1))
        return energy, res

    def fold(self, rna: Rna, cfg: EnergyCfg) -> tuple[Rna, CmdResult]:
        args = self._energy_cfg_args(cfg)
        with tempfile.NamedTemporaryFile("w") as f:
            if rna.r is None:
                raise ValueError(f"RNA {rna.name} has no sequence")
            f.write(rna.r)
            f.flush()
            res = self._run_cmd("./src/bin/RNAfold", *args, "--noPS", "-i", f.name)
            seq, db = res.stdout.strip().split("\n")
            db = db.split(" ")[0]
            predicted = RnaParser.parse(name=rna.name, seq=seq.strip(), db=db.strip())
        return predicted, res

    def partition(self, rna: Rna, cfg: EnergyCfg) -> None:
        raise NotImplementedError

    def subopt(
        self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg
    ) -> tuple[list[Rna], CmdResult]:
        args = self._energy_cfg_args(energy_cfg)
        args += self._subopt_cfg_args(subopt_cfg)
        res = self._run_cmd("./src/bin/RNAsubopt", *args, inp=rna.r)
        subopts = []
        for i in res.stdout.splitlines()[1:]:
            db, energy_str = re.split(r"\s+", i.strip())
            energy = Decimal(energy_str)
            subopts.append(Rna(name=rna.name, r=rna.r, s=db_to_secondary(db), energy=energy))
        return subopts, res

    def __str__(self) -> str:
        return "ViennaRNA"
