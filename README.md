# memerna

In all cases where an ordering of base_t p is used (e.g. data tables), it will be ACGU.

### Building

Run git submodule init and git submodule update to pull in external dependencies. 
Memerna requires a modern C++ compiler that supports C++14.

./build.py -t [debug, asan, ubsan, release, relwithdebinfo] -c

Will use clang with -c. The flag -t specifies the build type. Then run from ./build/memerna.

No guarantees this runs or even builds on Windows.

### Directories

- build: build output directory
- cmake: CMake scripts
- data: energy model data for memerna
- examples: various dot-bracket example folded RNAs
- extern: external projects and data (gtest, original data from rnastructure and nndb, rnark)
- scripts: scripts for various things (see below)
- src: source
- tests: tests

### Running the tests
Run from ./build/run_tests after building.

To run AFL, first build the afl binary with build.py -t relwithdebinfo -a, then run:

afl-fuzz -x extern/afl/dict.dct -m 2000 -t 1000 -i ./extern/afl/testcases 
  -o ./build/afl_test ./build/afl-clang-fast++-relwithdebinfo/fuzz -afl
  
You can optionally specify -random to fuzz.

### Useful commands

./scripts/run_benchmarks.py
Runs benchmarking for various packages.

./scripts/run.py --rnastructure-loc=<> -p <ct or db filename>
Supports EFN or folding through -e and -f respectively. Benchmarking through -b.
Can also specify memevault RNA via -kv <memevault name>.

./scripts/plot_benchmarks.py
Plots the output of the benchmarking.

./scripts/parse_data.py
Parses the original data files from orig_data and outputs them in memerna format in data/.

./scripts/extract_subsequence.py
Gives you the subsequence (0 indexed) between st, en of an RNA.

./scripts/crop.py
Crops images in a directory.

./scripts/convert.py
Converts between db and ct files.
