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


@dataclass
class ViennaRna(RnaPackage):
    def energy_cfg_args(self, cfg: EnergyCfg):
        args = []
        if not cfg.lonely_pairs:
            args.append("--no-lonely-pairs")
        match cfg.ctd:
            case CtdCfg.NONE:
                args.append("-d0")
            case CtdCfg.D2:
                args.append("-d2")
            case CtdCfg.NO_COAX:
                raise NotImplementedError(
                    "ViennaRNA does not support CTDs with no coaxial stacking",
                )
            case CtdCfg.CTD:
                args.append("-d3")
        return args

    def subopt_cfg_args(self, cfg: SuboptCfg):
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

    def efn(self, rna: Rna, cfg: EnergyCfg):
        args = self.energy_cfg_args(cfg)
        res = self._run_cmd("./src/bin/RNAeval", *args, input=f"{rna.r}\n{rna.db()}")
        match = re.search(r"\s+\(\s*([0-9\.\-]+)\s*\)", res.stdout.strip())
        energy = float(match.group(1))
        return energy, res

    def fold(self, rna: Rna, cfg: EnergyCfg):
        args = self.energy_cfg_args(cfg)
        with tempfile.NamedTemporaryFile("w") as f:
            f.write(rna.r)
            f.flush()
            res = self._run_cmd("./src/bin/RNAfold", *args, "--noPS", "-i", f.name)
            seq, db = res.stdout.strip().split("\n")
            db = db.split(" ")[0]
            predicted = Rna.from_name_seq_db(rna.name, seq.strip(), db.strip())
        return predicted, res

    def partition(self, rna: Rna, cfg: EnergyCfg):
        # TODO: Implement this.
        raise NotImplementedError

    def subopt(self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg):
        args = self.energy_cfg_args(energy_cfg)
        args += self.subopt_cfg_args(subopt_cfg)
        with tempfile.NamedTemporaryFile("r") as out:
            res = self._run_cmd("./src/bin/RNAsubopt", *args, input=rna.r)
            subopts = []
            for i in res.stdout.splitlines()[1:]:
                db, energy = re.split(r"\s+", i.strip())
                subopts.append(Rna(name=rna.name, r=rna.r, s=db_to_secondary(db), energy=energy))
        return subopts, res

    def __str__(self):
        return "ViennaRNA"
