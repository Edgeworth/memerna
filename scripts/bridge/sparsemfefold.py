# Copyright 2022 Eliot Courtney.


from dataclasses import dataclass
from pathlib import Path

from scripts.bridge.rnapackage import RnaPackage
from scripts.util.command import CmdLimits


@dataclass
class SparseMfeFold(RnaPackage):
    def efn(self, rna: Rna, cfg: EnergyCfg):
        raise NotImplementedError

    def fold(self, rna: Rna, cfg: EnergyCfg):
        res = run_command(
            os.path.join(self.path, "src", "SparseMFEFold"),
            input=rna.r,
            record_stdout=True,
        )
        seq, db = res.stdout.strip().split("\n")
        db = db.split(" ")[0]
        predicted = Rna.from_name_seq_db(rna.name, seq.strip(), db.strip())
        return predicted, res

    def partition(self, rna: Rna, cfg: EnergyCfg):
        raise NotImplementedError

    def subopt(self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg):
        raise NotImplementedError

    def __str__(self):
        return "SparseMFEFold"
