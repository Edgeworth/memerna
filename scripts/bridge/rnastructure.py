# Copyright 2022 E.


class RNAstructureDistribution:
    def __init__(self, loc=None):
        try:
            import scripts.default_paths

            loc = loc or scripts.default_paths.RNASTRUCTURE_PATH
        except ImportError:
            pass
        assert loc
        self.loc = fix_path(loc)
        os.putenv("DATAPATH", os.path.join(self.loc, "data_tables"))

    def fold(self, rna):
        with tempfile.NamedTemporaryFile("w") as f, tempfile.NamedTemporaryFile("r") as out:
            f.write(rna.to_seq_file())
            f.flush()
            res = run_command(os.path.join(self.loc, "exe", "Fold"), "-mfe", f.name, out.name)
            output = out.read()
            predicted = RNA.from_any_file(output)
        return predicted, res

    def efn(self, rna):
        with tempfile.NamedTemporaryFile("w") as f, tempfile.NamedTemporaryFile("r") as out:
            f.write(rna.to_ct_file())
            f.flush()
            # Note that not giving the -s flag doesn't make it logarithmic.
            # RNAstructure 5.8 adds the logarithmic and asymmetry models together in this case.
            # RNAstructure also uses a coefficient of -6 for the number of branches, rather than
            # the fitted -9.
            res = run_command(os.path.join(self.loc, "exe", "efn2"), "-s", f.name, out.name)
            output = out.read()
            match = re.search(r"[eE]nergy = (.+)", output.strip())
            energy = float(match.group(1))
        return energy, res

    def suboptimal(self, rna, delta, limits, num_only=False):
        with tempfile.NamedTemporaryFile("w") as f, tempfile.NamedTemporaryFile("r") as out:
            f.write(rna.to_seq_file())
            f.flush()
            res = try_command(
                os.path.join(self.loc, "exe", "AllSub"),
                "-a",
                f"{delta / 10.0:.2f}",
                f.name,
                out.name,
                limits=limits,
            )
            if num_only:
                res2 = run_command("wc", "-l", out.name, record_stdout=True)
                num_lines = int(res2.stdout.strip().split(" ")[0])
                retval = num_lines // (len(rna.seq) + 1)
            else:
                output = out.read()
                # TODO does not extract energy yet
                retval = [(0.0, i) for i in rnas_from_multi_ct_file(output)]
        return retval, res

    def close(self):
        pass

    def __str__(self):
        return "RNAstructureDistribution"


class RNAstructureHarness(HarnessFolder):
    def __init__(self, loc=None):
        super().__init__(loc, "RNAstructureHarness", "-r")
