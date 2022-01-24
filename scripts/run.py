#!/usr/bin/env python3
# Copyright 2016 Eliot Courtney.
import argparse
import shutil
import tempfile
import os

from scripts.common import *
from scripts.memevault import MemeVault
from scripts.rna import *


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
                "%.2f" % (delta / 10.0),
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
            os.path.join("harness"), "-f", self.flag, record_stdout=True, input=rna.seq
        )
        os.chdir(prev_dir)
        lines = res.stdout.strip().split("\n")
        assert len(lines) == 2
        return RNA.from_name_seq_db(rna.name, rna.seq, lines[1]), res

    def batch_efn(self, rnas):
        prev_dir = os.getcwd()
        os.chdir(self.loc)
        input = "\n".join("%s\n%s" % (rna.seq, rna.db()) for rna in rnas)
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
                        (float(energy) / 10.0, RNA.from_name_seq_db(rna.name, rna.seq, db))
                    )
            os.chdir(prev_dir)
        return retval, res

    def close(self):
        pass

    def __str__(self):
        return self.name


class RNAstructureHarness(HarnessFolder):
    def __init__(self, loc=None):
        super().__init__(loc, "RNAstructureHarness", "-r")


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
                limits=limits
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
                        (float(energy) / 10.0, RNA.from_name_seq_db(rna.name, rna.seq, db))
                    )
        return retval, res

    def close(self):
        pass

    def __str__(self):
        return self.name


class ViennaRNA:
    def __init__(self, loc=None, d3=False, sorted=False):
        try:
            import scripts.default_paths

            loc = loc or scripts.default_paths.VIENNARNA_PATH
        except ImportError:
            pass
        assert loc
        self.loc = fix_path(loc)
        self.extra_args = ["-d2"]
        self.name = "ViennaRNA-d2"
        if d3:
            self.extra_args = ["-d3"]
            self.name = "ViennaRNA-d3"
        if sorted:
            self.extra_args.append("--sorted")
            self.name += "-sorted"

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
                record_stdout=True
            )
            seq, db = res.stdout.strip().split("\n")
            db = db.split(" ")[0]
            predicted = RNA.from_name_seq_db(rna.name, seq.strip(), db.strip())
        return predicted, res

    def efn(self, rna):
        res = run_command(
            os.path.join(self.loc, "src", "bin", "RNAeval"),
            *self.extra_args,
            input=rna.seq + "\n" + rna.db(),
            record_stdout=True
        )
        match = re.search(r"\s+\(\s*([0-9\.\-]+)\s*\)", res.stdout.strip())  # mfw this regex
        energy = float(match.group(1))
        return energy, res

    def suboptimal(self, rna, delta, limits, num_only=False):
        with tempfile.NamedTemporaryFile("r") as out:
            res = try_command(
                os.path.join(self.loc, "src", "bin", "RNAsubopt"),
                *self.extra_args,
                "-e",
                "%.1f" % (delta / 10.0),
                input=rna.seq,
                record_stdout=out.name,
                limits=limits
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
                    retval.append((float(energy), RNA.from_name_seq_db(rna.name, rna.seq, db)))
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

    def fold(self, rna):
        raise NotImplementedError

    def efn(self, rna):
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
                "%.1f" % (delta / 10.0),
                input=rna.seq,
                record_stdout=out.name,
                limits=limits
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
                    retval.append((float(energy), RNA.from_name_seq_db(rna.name, rna.seq, db)))
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

    def fold(self, rna):
        raise NotImplementedError

    def efn(self, rna):
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
                "%.1f" % (delta / 10.0),
                "-input",
                f.name,
                limits=limits
            )
            files = ["subopt-%d.stdout" % i for i in range(self.n)]
            if num_only:
                retval = 0
                for file in files:
                    res2 = run_command("wc", "-l", file, record_stdout=True)
                    num_lines = int(res2.stdout.strip().split(" ")[0])
                    retval += num_lines - 2
            else:
                retval = []
                for file in files:
                    output = open(file, "r").read()
                    for i in output.splitlines()[2:]:
                        db, energy = re.split(r"\s+", i.strip())[:2]
                        retval.append((float(energy), RNA.from_name_seq_db(rna.name, rna.seq, db)))
        os.chdir(prev_dir)
        return retval, res

    def close(self):
        shutil.rmtree(self.tempdir)

    def __str__(self):
        return self.name


class UNAFold:
    def __init__(self, loc=None):
        try:
            import scripts.default_paths

            loc = loc or scripts.default_paths.UNAFOLD_PATH
        except ImportError:
            pass
        assert loc
        self.loc = fix_path(loc)
        self.tempdir = tempfile.mkdtemp()
        os.putenv("UNAFOLDDAT", os.path.join(self.loc, "data"))

    def fold(self, rna):
        prev_dir = os.getcwd()
        os.chdir(self.tempdir)
        with open(os.path.join(self.tempdir, "rna.seq"), "w") as f:
            f.write(rna.to_db_file())
            f.flush()
            res = run_command(os.path.join(self.loc, "src", "hybrid-ss-min"), f.name)
            predicted = RNA.from_any_file(read_file(os.path.splitext(f.name)[0] + ".ct"))
        os.chdir(prev_dir)
        return predicted, res

    def efn(self, rna):
        prev_dir = os.getcwd()
        os.chdir(self.tempdir)
        with open(os.path.join(self.tempdir, "rna.seq"), "w") as f:
            f.write(rna.to_ct_file())
            f.flush()
            res = run_command(
                os.path.join(self.loc, "src", "ct-energy"), f.name, record_stdout=True
            )
            energy = float(res.stdout.strip())
        os.chdir(prev_dir)
        return energy, res

    def suboptimal(self, rna, delta, limits, num_only=False):
        raise NotImplementedError

    def close(self):
        shutil.rmtree(self.tempdir)

    def __str__(self):
        return "UNAFold"


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
            os.path.join(self.loc, "src", "SparseMFEFold"), input=rna.seq, record_stdout=True
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


def run_fold(program, rna):
    frna, res = program.fold(rna)
    print("Folding %s with %s: %s\n  %s" % (rna.name, program, frna.db(), res))


def run_efn(program, rna):
    energy, res = program.efn(rna)
    print("Energy of %s with %s: %f\n  %s" % (rna.name, program, energy, res))


def run_suboptimal(program, rna, delta, subopt_max_print):
    subopts, res = program.suboptimal(rna, delta, None)
    if res.ret:
        print("Execution failed")
        return
    print("%d suboptimal structures of %s with %s - %s" % (len(subopts), rna.name, program, res))
    subopt_subset = subopts
    if subopt_max_print > 0:
        subopt_subset = subopts[:subopt_max_print]
    for energy, structure in subopt_subset:
        print("%f %s" % (energy, structure.db()))


BENCHMARK_NUM_TRIES = 5
# 12 GiB, 4 minutes
BENCHMARK_LIMITS = (12 * 1024 * 1024, 5 * 60)


def run_subopt_benchmark(program, dataset, delta, num_file):
    nums = None
    if num_file:
        assert isinstance(program, MemeRNA)  # Only memerna supports num folding.
        nums = {a: int(b) for a, b in [i.split(" ") for i in read_file(num_file).splitlines()]}
    memevault = MemeVault(dataset)
    print("Benchmarking suboptimals with %s on %s with delta %d" % (program, dataset, delta))
    filename = "%s_%s_subopt_%d.results" % (program, dataset, delta)
    if os.path.exists(filename):
        print("Not overwriting %s" % filename)
        sys.exit(1)
    with open(filename, "w") as f:
        idx = 1
        for rna in memevault:
            print("Running %s on #%d %s" % (program, idx, rna.name))
            idx += 1
            len_res = []
            failed = False
            for i in range(BENCHMARK_NUM_TRIES):
                if nums:
                    length, res = program.suboptimal_maxdelta(
                        rna, nums[rna.name], -1, BENCHMARK_LIMITS, True
                    )
                else:
                    length, res = program.suboptimal(rna, delta, BENCHMARK_LIMITS, True)
                if res.ret:
                    print("Got OOM or OOT (presumably), skipping.")
                    failed = True
                    break
                len_res.append((length, res))
            if failed:
                continue
            for i, lr in enumerate(len_res):
                num_subopt, res = lr
                f.write(
                    "%s %d %d %.5f %.5f %.5f %d\n"
                    % (rna.name, i, len(rna.seq), res.real, res.usersys, res.maxrss, num_subopt)
                )


def run_fold_benchmark(program, dataset, rnastructure_harness):
    if rnastructure_harness is None:
        print("Need RNAstructure for benchmark efn checking")
        sys.exit(1)
    memevault = MemeVault(dataset)
    print("Benchmarking folding with %s on %s" % (program, dataset))
    filename = "%s_%s_fold.results" % (program, dataset)
    if os.path.exists(filename):
        print("Not overwriting %s" % filename)
        sys.exit(1)
    with open(filename, "w") as f:
        idx = 1
        for rna in memevault:
            print("Running %s on #%d %s" % (program, idx, rna.name))
            prs = []
            for i in range(BENCHMARK_NUM_TRIES):
                prs.append(program.fold(rna))

            accuracy = RNAAccuracy.from_rna(rna, prs[0][0])
            energy, _ = rnastructure_harness.efn(prs[0][0])
            for i, pr in enumerate(prs):
                predicted, res = pr
                f.write(
                    "%s %d %d %.5f %.5f %.5f %.5f %.5f %.5f %.2f\n"
                    % (
                        rna.name,
                        i,
                        len(rna.seq),
                        res.real,
                        res.usersys,
                        res.maxrss,
                        accuracy.fscore,
                        accuracy.ppv,
                        accuracy.sensitivity,
                        energy,
                    )
                )
            idx += 1


def parse_rna_from_args(parser, args):
    if args.benchmark:
        return None
    if bool(args.path) + bool(args.memevault) + bool(args.cmd) != 1:
        parser.error("Exactly one primary/secondary sequence required")

    if args.path:
        return RNA.from_any_file(read_file(args.path))
    elif args.memevault:
        return MemeVault("archiveii")[args.memevault]
    elif args.cmd:
        if args.fold or args.subopt:
            if len(args.cmd) != 1:
                parser.error("Direct specification requires one argument for prediction.")
            seq = args.cmd[0]
            return RNA("cmdline", seq, [-1] * len(seq))
        elif args.energy:
            if len(args.cmd) != 2:
                parser.error("Direct specification requires two arguments for efn.")
            return RNA.from_name_seq_db("cmdline", *args.cmd)


def process_command(*extra_args):
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", type=str)
    parser.add_argument("-mv", "--memevault", type=str)
    parser.add_argument("cmd", type=str, nargs="*")

    parser.add_argument("--rnastructure-loc")
    parser.add_argument("--memerna-loc")
    parser.add_argument("--viennarna-loc")
    parser.add_argument("--unafold-loc")
    parser.add_argument("--sparsemfefold-loc")
    parser.add_argument("--sjsvienna-loc")
    parser.add_argument("--sjsviennampi-loc")

    parser.add_argument("-rh", "--rnastructure-harness", action="store_true")
    parser.add_argument("-rd", "--rnastructure-distribution", action="store_true")
    parser.add_argument("-vd2", "--viennarna-d2", action="store_true")
    parser.add_argument("-vd3", "--viennarna-d3", action="store_true")
    parser.add_argument("-vd2s", "--viennarna-d2-sorted", action="store_true")
    parser.add_argument("-vd3s", "--viennarna-d3-sorted", action="store_true")
    parser.add_argument("-sjs", "--sjsvienna", action="store_true")
    parser.add_argument("-sjss", "--sjsvienna-sorted", action="store_true")
    parser.add_argument("-sjsmpi", "--sjsviennampi", action="store_true")
    parser.add_argument("-sjsmpis", "--sjsviennampi-sorted", action="store_true")
    parser.add_argument("-u", "--unafold", action="store_true")
    parser.add_argument("-m", "--memerna", action="store_true")
    parser.add_argument("-ms", "--memerna-sorted", action="store_true")
    parser.add_argument("-smf", "--sparsemfefold", action="store_true")

    parser.add_argument("-f", "--fold", action="store_true")
    parser.add_argument("-e", "--energy", action="store_true")
    parser.add_argument("-s", "--subopt", type=int)
    parser.add_argument("-b", "--benchmark", type=str)
    parser.add_argument("-bn", "--benchmark_nums", type=str)
    parser.add_argument("--subopt-max-print", type=int, default=20)

    args = parser.parse_args(sys.argv[1:] + list(*extra_args))

    if bool(args.fold) + bool(args.energy) + bool(args.subopt is not None) != 1:
        parser.error("Exactly one of --fold, --energy, or --subopt is required.")

    programs = []
    rnastructure_harness = RNAstructureHarness(args.memerna_loc)
    if args.rnastructure_harness:
        programs.append(rnastructure_harness)
    if args.rnastructure_distribution:
        programs.append(RNAstructureDistribution(args.rnastructure_loc))
    if args.viennarna_d2:
        programs.append(ViennaRNA(args.viennarna_loc, False))
    if args.viennarna_d3:
        programs.append(ViennaRNA(args.viennarna_loc, True))
    if args.viennarna_d2_sorted:
        programs.append(ViennaRNA(args.viennarna_loc, False, True))
    if args.viennarna_d3_sorted:
        programs.append(ViennaRNA(args.viennarna_loc, True, True))
    if args.sjsvienna:
        programs.append(SJSVienna(args.sjsvienna_loc, False))
    if args.sjsvienna_sorted:
        programs.append(SJSVienna(args.sjsvienna_loc, True))
    if args.sjsviennampi:
        programs.append(SJSViennaMPI(args.sjsviennampi_loc, False))
    if args.sjsviennampi_sorted:
        programs.append(SJSViennaMPI(args.sjsviennampi_loc, True))
    if args.unafold:
        programs.append(UNAFold(args.unafold_loc))
    if args.memerna:
        programs.append(MemeRNA(args.memerna_loc, False))
    if args.memerna_sorted:
        programs.append(MemeRNA(args.memerna_loc, True))
    if args.sparsemfefold:
        programs.append(SparseMFEFold(args.sparsemfefold_loc))

    rna = parse_rna_from_args(parser, args)
    for program in programs:
        if args.benchmark:
            if args.fold:
                run_fold_benchmark(program, args.benchmark, rnastructure_harness)
            elif args.subopt:
                run_subopt_benchmark(program, args.benchmark, args.subopt, args.benchmark_nums)
        elif args.fold:
            run_fold(program, rna)
        elif args.energy:
            run_efn(program, rna)
        elif args.subopt:
            run_suboptimal(program, rna, args.subopt, args.subopt_max_print)

    for program in programs:
        program.close()


if __name__ == "__main__":
    process_command()
