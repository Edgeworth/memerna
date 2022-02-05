# Copyright 2022 Eliot Courtney.


from pathlib import Path


class ViennaRna(RnaPackage):
    def __init__(self, loc: Path, d3=False, sorted=False):
        self.loc = loc
        self.extra_args = ["-d2"]
        self.name = "ViennaRNA-d2"
        if d3:
            self.extra_args = ["-d3"]
            self.name = "ViennaRNA-d3"
        if sorted:
            self.extra_args.append("--sorted")
            self.name += "-sorted"

    def efn(self, rna):
        res = run_command(
            os.path.join(self.loc, "src", "bin", "RNAeval"),
            *self.extra_args,
            input=f"{rna.seq}\n{rna.db()}",
            record_stdout=True,
        )
        match = re.search(r"\s+\(\s*([0-9\.\-]+)\s*\)", res.stdout.strip())  # mfw this regex
        energy = float(match.group(1))
        return energy, res

    def fold(self, rna):
        with tempfile.NamedTemporaryFile("w") as f:
            f.write(rna.seq)
            f.flush()
            res = run_command(
                os.path.join(self.loc, "src", "bin", "RNAfold"),
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

    def suboptimal(self, rna, delta, limits, num_only=False):
        with tempfile.NamedTemporaryFile("r") as out:
            res = try_command(
                os.path.join(self.loc, "src", "bin", "RNAsubopt"),
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


class SJSVienna:
    def __init__(self, loc=None, sorted=False):
        try:
            import scripts.default_paths

            loc = loc or scripts.default_paths.SJSVIENNA_PATH
        except ImportError:
            pass
        assert loc
        self.loc = fix_path(loc)
        self.extra_args = []
        self.name = "SJSVienna"
        if sorted:
            self.extra_args.append("-s")
            self.name += "-sorted"

    def efn(self, rna):
        raise NotImplementedError

    def fold(self, rna):
        raise NotImplementedError

    def suboptimal(self, rna, delta, limits, num_only=False):
        with tempfile.NamedTemporaryFile("r") as out:
            res = try_command(
                os.path.join(self.loc, "Progs", "RNAsubopt"),
                "-P",
                os.path.join(self.loc, "2004.par"),
                *self.extra_args,
                "-d2",
                "-e",
                f"{delta / 10.0:.1f}",
                input=rna.seq,
                record_stdout=out.name,
                limits=limits,
            )
            if num_only:
                res2 = run_command("wc", "-l", out.name, record_stdout=True)
                num_lines = int(res2.stdout.strip().split(" ")[0])
                retval = num_lines - 2
            else:
                retval = []
                output = out.read()
                for i in output.splitlines()[2:]:
                    db, energy = re.split(r"\s+", i.strip())[:2]
                    retval.append((float(energy), Rna.from_name_seq_db(rna.name, rna.seq, db)))
        return retval, res

    def close(self):
        pass

    def __str__(self):
        return self.name


class SJSViennaMPI:
    def __init__(self, loc=None, sorted=False, n=4):
        try:
            import scripts.default_paths

            loc = loc or scripts.default_paths.SJSVIENNAMPI_PATH
        except ImportError:
            pass
        assert loc
        self.loc = fix_path(loc)
        self.extra_args = []
        self.name = "SJSViennaMPI"
        self.n = n
        self.tempdir = tempfile.mkdtemp()
        if sorted:
            self.extra_args.append("-s")
            self.name += "-sorted"

    def efn(self, rna):
        raise NotImplementedError

    def fold(self, rna):
        raise NotImplementedError

    def suboptimal(self, rna, delta, limits, num_only=False):
        prev_dir = os.getcwd()
        os.chdir(self.tempdir)
        with tempfile.NamedTemporaryFile("w") as f:
            f.write(rna.seq)
            f.flush()
            res = try_command(
                "mpirun",
                "-n",
                str(self.n),
                os.path.join(self.loc, "Progs", "RNAsubopt"),
                "-P",
                os.path.join(self.loc, "2004.par"),
                *self.extra_args,
                "-d2",
                "-e",
                f"{delta / 10.0:.1f}",
                "-input",
                f.name,
                limits=limits,
            )
            files = [f"subopt-{int(i)}.stdout" for i in range(self.n)]
            if num_only:
                retval = 0
                for file in files:
                    res2 = run_command("wc", "-l", file, record_stdout=True)
                    num_lines = int(res2.stdout.strip().split(" ")[0])
                    retval += num_lines - 2
            else:
                retval = []
                for file in files:
                    output = open(file).read()
                    for i in output.splitlines()[2:]:
                        db, energy = re.split(r"\s+", i.strip())[:2]
                        retval.append((float(energy), Rna.from_name_seq_db(rna.name, rna.seq, db)))
        os.chdir(prev_dir)
        return retval, res

    def close(self):
        shutil.rmtree(self.tempdir)

    def __str__(self):
        return self.name
