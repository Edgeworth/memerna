#!/usr/bin/env python3
# Copyright 2022 E.
import argparse
import os
import shutil
import tempfile

from scripts.common import *
from scripts.memevault import MemeVault
from scripts.model.rna import *


def run_fold(program, rna):
    frna, res = program.fold(rna)
    print(f"Folding {rna.name} with {program}: {frna.db()}\n  {res}")


def run_efn(program, rna):
    energy, res = program.efn(rna)
    print(f"Energy of {rna.name} with {program}: {energy:f}\n  {res}")


def run_suboptimal(program, rna, delta, subopt_max_print):
    subopts, res = program.suboptimal(rna, delta, None)
    if res.ret:
        print("Execution failed")
        return
    print(f"{len(subopts)} suboptimal structures of {rna.name} with {program} - {res}")
    subopt_subset = subopts
    if subopt_max_print > 0:
        subopt_subset = subopts[:subopt_max_print]
    for energy, structure in subopt_subset:
        print(f"{energy:f} {structure.db()}")


BENCHMARK_NUM_TRIES = 5
# 12 GiB, 4 minutes
BENCHMARK_LIMITS = (12 * 1024 * 1024, 5 * 60)


def run_subopt_benchmark(program, dataset, delta, num_file):
    nums = None
    if num_file:
        assert isinstance(program, MemeRNA)  # Only memerna supports num folding.
        nums = {a: int(b) for a, b in [i.split(" ") for i in read_file(num_file).splitlines()]}
    memevault = MemeVault(dataset)
    print(f"Benchmarking suboptimals with {program} on {dataset} with delta {int(delta)}")
    filename = f"{program}_{dataset}_subopt_{int(delta)}.results"
    if os.path.exists(filename):
        print(f"Not overwriting {filename}")
        sys.exit(1)
    with open(filename, "w") as f:
        idx = 1
        for rna in memevault:
            print(f"Running {program} on #{int(idx)} {rna.name}")
            idx += 1
            len_res = []
            failed = False
            for i in range(BENCHMARK_NUM_TRIES):
                if nums:
                    length, res = program.suboptimal_maxdelta(
                        rna,
                        nums[rna.name],
                        -1,
                        BENCHMARK_LIMITS,
                        True,
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
                    f"{rna.name} {int(i)} {len(rna.seq)} {res.real:.5f} {res.usersys:.5f} {res.maxrss:.5f} {int(num_subopt)}\n",
                )


def run_fold_benchmark(program, dataset, rnastructure_harness):
    if rnastructure_harness is None:
        print("Need RNAstructure for benchmark efn checking")
        sys.exit(1)
    memevault = MemeVault(dataset)
    print(f"Benchmarking folding with {program} on {dataset}")
    filename = f"{program}_{dataset}_fold.results"
    if os.path.exists(filename):
        print(f"Not overwriting {filename}")
        sys.exit(1)
    with open(filename, "w") as f:
        idx = 1
        for rna in memevault:
            print(f"Running {program} on #{int(idx)} {rna.name}")
            prs = []
            for i in range(BENCHMARK_NUM_TRIES):
                prs.append(program.fold(rna))

            accuracy = RNAAccuracy.from_rna(rna, prs[0][0])
            energy, _ = rnastructure_harness.efn(prs[0][0])
            for i, pr in enumerate(prs):
                predicted, res = pr
                f.write(
                    f"{rna.name} {int(i)} {len(rna.seq)} {res.real:.5f} {res.usersys:.5f} {res.maxrss:.5f} {accuracy.fscore:.5f} {accuracy.ppv:.5f} {accuracy.sensitivity:.5f} {energy:.2f}\n",
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
