# Copyright 2022 Eliot Courtney.


from pathlib import Path


class ViennaRna(RnaPackage):
    def __init__(self, path: Path):
        self.path = path

    def efn(self, rna: Rna, cfg: EnergyCfg):
        res = run_command(
            os.path.join(self.path, "src", "bin", "RNAeval"),
            *self.extra_args,
            input=f"{rna.seq}\n{rna.db()}",
            record_stdout=True,
        )
        match = re.search(r"\s+\(\s*([0-9\.\-]+)\s*\)", res.stdout.strip())  # mfw this regex
        energy = float(match.group(1))
        return energy, res

    def fold(self, rna: Rna, cfg: EnergyCfg):
        with tempfile.NamedTemporaryFile("w") as f:
            f.write(rna.seq)
            f.flush()
            res = run_command(
                os.path.join(self.path, "src", "bin", "RNAfold"),
                *self.extra_args,
                "--noPS",
                "-i",
                f.name,
                record_stdout=True,
            )
            seq, db = res.stdout.strip().split("\n")
            db = db.split(" ")[0]
            predicted = Rna.from_name_seq_db(rna.name, seq.strip(), db.strip())
        return predicted, res

    def partition(self, rna: Rna, cfg: EnergyCfg):
        raise NotImplementedError

    def subopt(self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg):
        with tempfile.NamedTemporaryFile("r") as out:
            res = try_command(
                os.path.join(self.path, "src", "bin", "RNAsubopt"),
                *self.extra_args,
                "-e",
                f"{delta / 10.0:.1f}",
                input=rna.seq,
                record_stdout=out.name,
                limits=limits,
            )
            if num_only:
                res2 = run_command("wc", "-l", out.name, record_stdout=True)
                num_lines = int(res2.stdout.strip().split(" ")[0])
                retval = num_lines - 1
            else:
                retval = []
                output = out.read()
                for i in output.splitlines()[1:]:
                    db, energy = re.split(r"\s+", i.strip())
                    retval.append((float(energy), Rna.from_name_seq_db(rna.name, rna.seq, db)))
        return retval, res

    def close(self):
        pass

    def __str__(self):
        return self.name
