# Copyright 2022 Eliot Courtney.


class SparseMFEFold:
    def __init__(self, loc=None):
        try:
            import scripts.default_paths

            loc = loc or scripts.default_paths.SPARSEMFEFOLD_PATH
        except ImportError:
            pass
        assert loc
        self.loc = fix_path(loc)

    def fold(self, rna):
        res = run_command(
            os.path.join(self.loc, "src", "SparseMFEFold"),
            input=rna.seq,
            record_stdout=True,
        )
        seq, db = res.stdout.strip().split("\n")
        db = db.split(" ")[0]
        predicted = RNA.from_name_seq_db(rna.name, seq.strip(), db.strip())
        return predicted, res

    def efn(self, rna):
        raise NotImplementedError

    def suboptimal(self, rna, delta, limits, num_only=False):
        raise NotImplementedError

    def close(self):
        pass

    def __str__(self):
        return "SparseMFEFold"
