# Copyright 2022 E.


from pathlib import Path


class SparseMfeFold(RnaPackage):
    def __init__(self, loc: Path):
        self.loc = loc

    def efn(self, rna: Rna, cfg: EnergyCfg):
        raise NotImplementedError

    def fold(self, rna: Rna, cfg: EnergyCfg):
        res = run_command(
            os.path.join(self.loc, "src", "SparseMFEFold"),
            input=rna.seq,
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
