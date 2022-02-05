# Copyright 2022 Eliot Courtney.


from pathlib import Path


# TODO: create tempdir separately and delete it after each use
class UnaFold(RnaPackage):
    def __init__(self, loc: Path):
        self.loc = loc
        self.tempdir = tempfile.mkdtemp()
        os.putenv("UNAFOLDDAT", os.path.join(self.loc, "data"))

    def efn(self, rna):
        prev_dir = os.getcwd()
        os.chdir(self.tempdir)
        with open(os.path.join(self.tempdir, "rna.seq"), "w") as f:
            f.write(rna.to_ct_file())
            f.flush()
            res = run_command(
                os.path.join(self.loc, "src", "ct-energy"),
                f.name,
                record_stdout=True,
            )
            energy = float(res.stdout.strip())
        os.chdir(prev_dir)
        return energy, res

    def fold(self, rna):
        prev_dir = os.getcwd()
        os.chdir(self.tempdir)
        with open(os.path.join(self.tempdir, "rna.seq"), "w") as f:
            f.write(rna.to_db_file())
            f.flush()
            res = run_command(os.path.join(self.loc, "src", "hybrid-ss-min"), f.name)
            predicted = Rna.from_any_file(read_file(f"{os.path.splitext(f.name)[0]}.ct"))
        os.chdir(prev_dir)
        return predicted, res

    def suboptimal(self, rna, delta, limits, num_only=False):
        raise NotImplementedError

    def close(self):
        shutil.rmtree(self.tempdir)

    def __str__(self):
        return "UNAFold"
