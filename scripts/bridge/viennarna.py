# Copyright 2022 Eliot Courtney.
from dataclasses import dataclass
import re
import tempfile

from scripts.bridge.rnapackage import RnaPackage
from scripts.model.config import CtdCfg
from scripts.model.config import EnergyCfg
from scripts.model.config import SuboptCfg
from scripts.model.parse import db_to_secondary
from scripts.model.rna import Rna
from scripts.util.command import CmdResult
from typing import Tuple


@dataclass
class ViennaRna(RnaPackage):
    def energy_cfg_args(self, cfg: EnergyCfg) -> list[str]:
        args = []
        if not cfg.lonely_pairs:
            args.append("--no-lonely-pairs")
        if cfg.ctd == CtdCfg.NONE:
            args.append("-d0")
        if cfg.ctd == CtdCfg.D2:
            args.append("-d2")
        if cfg.ctd == CtdCfg.NO_COAX:
            raise NotImplementedError(
                "ViennaRNA does not support CTDs with no coaxial stacking",
            )
        if cfg.ctd == CtdCfg.ALL:
            args.append("-d3")
        return args

    def subopt_cfg_args(self, cfg: SuboptCfg) -> list[str]:
        args = []
        if cfg.delta:
            args += ["--deltaEnergy", f"{cfg.delta / 10.0:.1f}"]
        if cfg.strucs:
            raise NotImplementedError(
                "ViennaRNA does not support reporting a maximum number of suboptimal structures",
            )
        if cfg.sorted:
            args += ["--sorted"]
        return args

    def efn(self, rna: Rna, cfg: EnergyCfg) -> Tuple[float, CmdResult]:
        args = self.energy_cfg_args(cfg)
        res = self._run_cmd("./src/bin/RNAeval", *args, inp=f"{rna.r}\n{rna.db()}")
        match = re.search(r"\s+\(\s*([0-9\.\-]+)\s*\)", res.stdout.strip())
        assert match is not None
        energy = float(match.group(1))
        return energy, res

    def fold(self, rna: Rna, cfg: EnergyCfg) -> Tuple[Rna, CmdResult]:
        args = self.energy_cfg_args(cfg)
        with tempfile.NamedTemporaryFile("w") as f:
            assert rna.r is not None
            f.write(rna.r)
            f.flush()
            res = self._run_cmd("./src/bin/RNAfold", *args, "--noPS", "-i", f.name)
            seq, db = res.stdout.strip().split("\n")
            db = db.split(" ")[0]
            predicted = Rna.parse(name=rna.name, seq=seq.strip(), db=db.strip())
        return predicted, res

    def partition(self, rna: Rna, cfg: EnergyCfg) -> None:
        raise NotImplementedError

    def subopt(
        self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg
    ) -> Tuple[list[Rna], CmdResult]:
        args = self.energy_cfg_args(energy_cfg)
        args += self.subopt_cfg_args(subopt_cfg)
        res = self._run_cmd("./src/bin/RNAsubopt", *args, inp=rna.r)
        subopts = []
        for i in res.stdout.splitlines()[1:]:
            db, energy_str = re.split(r"\s+", i.strip())
            energy = int(float(energy_str) * 10)
            subopts.append(Rna(name=rna.name, r=rna.r, s=db_to_secondary(db), energy=energy))
        return subopts, res

    def __str__(self) -> str:
        return "ViennaRNA"
