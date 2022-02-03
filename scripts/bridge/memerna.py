# Copyright 2022 Eliot Courtney.


class MemeRNA:
    def __init__(self, loc=None, sorted=False, ctd_output=False):
        try:
            import scripts.default_paths

            loc = loc or scripts.default_paths.MEMERNA_PATH
        except ImportError:
            pass
        assert loc
        self.loc = fix_path(loc)
        self.subopt_args = []
        self.name = "memerna"
        self.ctd_output = ctd_output
        if sorted:
            self.subopt_args.append("-sorted")
            self.name += "-sorted"
        if ctd_output:
            self.subopt_args = ["-ctd-output"]
            self.name += "-ctd"

    def fold(self, rna):
        prev_dir = os.getcwd()
        os.chdir(self.loc)
        res = run_command(os.path.join("fold"), rna.seq, record_stdout=True)
        os.chdir(prev_dir)
        _, db, _ = res.stdout.strip().split("\n")
        predicted = RNA.from_name_seq_db(rna.name, rna.seq, db.strip())
        return predicted, res

    def efn(self, rna):
        prev_dir = os.getcwd()
        os.chdir(self.loc)
        res = run_command(os.path.join("efn"), rna.seq, rna.db(), record_stdout=True)
        os.chdir(prev_dir)
        energy = float(res.stdout.strip().split(" ")[1]) / 10.0
        return energy, res

    def suboptimal(self, rna, delta, limits, num_only=False):
        return self.suboptimal_maxdelta(rna, -1, delta, limits, num_only)

    def suboptimal_maxdelta(self, rna, max_num, delta, limits, num_only):
        maxdelta_args = []
        if delta >= 0:
            maxdelta_args += ["-delta", str(delta)]
        if max_num >= 0:
            maxdelta_args += ["-num", str(max_num)]
        assert len(maxdelta_args) > 0
        with tempfile.NamedTemporaryFile("r") as out:
            prev_dir = os.getcwd()
            os.chdir(self.loc)
            res = try_command(
                os.path.join("subopt"),
                *self.subopt_args,
                *maxdelta_args,
                rna.seq,
                record_stdout=out.name,
                limits=limits,
            )
            os.chdir(prev_dir)
            if num_only:
                res2 = run_command("wc", "-l", out.name, record_stdout=True)
                num_lines = int(res2.stdout.strip().split(" ")[0])
                retval = num_lines - 1
            else:
                assert not self.ctd_output
                retval = []
                output = out.read()
                for i in output.splitlines()[:-1]:
                    energy, db = re.split(r"\s+", i.strip())
                    retval.append(
                        (float(energy) / 10.0, RNA.from_name_seq_db(rna.name, rna.seq, db)),
                    )
        return retval, res

    def close(self):
        pass

    def __str__(self):
        return self.name


class HarnessFolder:
    def __init__(self, loc, name, flag):
        try:
            import scripts.default_paths

            loc = loc or scripts.default_paths.MEMERNA_PATH
        except ImportError:
            pass
        assert loc
        self.loc = fix_path(loc)
        self.name = name
        self.flag = flag

    def fold(self, rna):
        prev_dir = os.getcwd()
        os.chdir(self.loc)
        res = run_command(
            os.path.join("harness"),
            "-f",
            self.flag,
            record_stdout=True,
            input=rna.seq,
        )
        os.chdir(prev_dir)
        lines = res.stdout.strip().split("\n")
        assert len(lines) == 2
        return RNA.from_name_seq_db(rna.name, rna.seq, lines[1]), res

    def batch_efn(self, rnas):
        prev_dir = os.getcwd()
        os.chdir(self.loc)
        input = "\n".join(f"{rna.seq}\n{rna.db()}" for rna in rnas)
        res = run_command(os.path.join("harness"), "-e", self.flag, record_stdout=True, input=input)
        os.chdir(prev_dir)
        energies = [float(i) / 10.0 for i in res.stdout.strip().split("\n")]
        return energies, res

    def efn(self, rna):
        energies, res = self.batch_efn([rna])
        return energies[0], res

    def suboptimal(self, rna, delta, limits, num_only=False):
        with tempfile.NamedTemporaryFile("r") as out:
            prev_dir = os.getcwd()
            os.chdir(self.loc)
            res = try_command(
                os.path.join("harness"),
                "-f",
                self.flag,
                "-subopt-delta",
                str(delta),
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
                for i in output.splitlines()[:-1]:
                    energy, db = re.split(r"\s+", i.strip())
                    retval.append(
                        (float(energy) / 10.0, RNA.from_name_seq_db(rna.name, rna.seq, db)),
                    )
            os.chdir(prev_dir)
        return retval, res

    def close(self):
        pass

    def __str__(self):
        return self.name
