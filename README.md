# memerna
Check out docs/thesis/thesis.pdf for my thesis which will explain a bit about
what this is all about.

In all cases where an ordering of base_t p is used (e.g. data tables), it will be ACGU.

### Building
The rest of this document uses $MRNA to locate the memerna directory.

Run git submodule init and git submodule update to pull in external dependencies.
Memerna requires a modern C++ compiler that supports C++14.

./build.py -t [debug, asan, ubsan, release, relwithdebinfo] -c

Will use clang with -c. The flag -t specifies the build type. Then run from $PREFIX/memerna.

No guarantees this runs or even builds on Windows.

### Directories

- build: build output directory
- cmake: CMake scripts
- data: energy model data for memerna
- docs: documentation
- examples: various dot-bracket example folded RNAs
- extern: external projects and data (original data from rnastructure and nndb, rnark)
- scripts: scripts for various things (see below)
- src: source
- tests: tests

### Running the tests
Run from $MRNA/run_tests after building.

### Fuzzing

#### Randomized fuzzing
```
make -j32 && ./fuzz -rnastructure-data $MRNA/extern/miles_rnastructure/data_tables/ \
  -memerna-data $MRNA/data/ 6 8 -print-interval 5 -no-partition -no-subopt
```

Use the -no-table-check option to only compare the result of memerna vs another
program, rather than the internal dp tables.

#### AFL
To run AFL++, first build the afl binary with build.py -t relwithdebinfo -a, then run:

```
sudo sh -c 'echo core >/proc/sys/kernel/core_pattern'
AFL_AUTORESUME=1 AFL_IMPORT_FIRST=1 AFL_TESTCACHE_SIZE=500 AFL_SKIP_CPUFREQ=1 \
  afl-fuzz -x $MRNA/extern/afl/fuzz/dict.dct -m 2000 -t 2000 \
  -i $MRNA/extern/afl/fuzz/testcases -o ./afl -- ./fuzz -afl \
  -memerna-data $MRNA/data/ -rnastructure-data $MRNA/extern/miles_rnastructure/data_tables/ \
  -rnastructure -table-check
```

Minimising test cases:
```
afl-tmin -i case -o ./afl/min -- ./fuzz -afl -memerna-data $MRNA/data/ \
  -rnastructure-data $MRNA/extern/miles_rnastructure/data_tables/ \
  -rnastructure -table-check
```

### Useful commands

./scripts/run_benchmarks.py
Runs benchmarking for various packages.

./scripts/run.py
Supports EFN or folding through -e and -f respectively. Benchmarking through -b.
Can also specify memevault RNA via -kv <memevault name>.

./scripts/parse_data.py
Parses the original data files from orig_data and outputs them in memerna format in data/.

./scripts/extract_subsequence.py
Gives you the subsequence (0 indexed) between st, en of an RNA.

./scripts/crop.py
Crops images in a directory.

./scripts/convert.py
Converts between db and ct files.
