# Copyright 2022 Eliot Courtney.
BENCHMARK_NUM_TRIES = 5
# 12 GiB, 4 minutes
BENCHMARK_LIMITS = (12 * 1024 * 1024 * 1024, 5 * 60)


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
                    f"{rna.name} {int(i)} {len(rna.r)} {res.real:.5f} {res.usersys:.5f} {res.maxrss:.5f} {int(num_subopt)}\n",
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
                    f"{rna.name} {int(i)} {len(rna.r)} {res.real:.5f} {res.usersys:.5f} {res.maxrss:.5f} {accuracy.fscore:.5f} {accuracy.ppv:.5f} {accuracy.sensitivity:.5f} {energy:.2f}\n",
                )
            idx += 1


def benchmark():
    init_package_limits(**fn_args())
