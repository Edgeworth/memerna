# Copyright 2022 Eliot Courtney.


from pathlib import Path


class SparseMfeFold(RnaPackage):
    def __init__(self, loc: Path):
        self.loc = loc

    def efn(self, rna):
        raise NotImplementedError

    def fold(self, rna):
        res = run_command(
            os.path.join(self.loc, "src", "SparseMFEFold"),
            input=rna.seq,
            record_stdout=True,
        )
        seq, db = res.stdout.strip().split("\n")
        db = db.split(" ")[0]
        predicted = Rna.from_name_seq_db(rna.name, seq.strip(), db.strip())
        return predicted, res

    def suboptimal(self, rna, delta, limits, num_only=False):
        raise NotImplementedError

    def __str__(self):
        return "SparseMFEFold"
