# memerna
Check out docs/thesis/thesis.pdf for my thesis which will explain a bit about
what this is all about.

In all cases where an ordering of base_t p is used (e.g. data tables), it will be ACGU.

### Building
The rest of this document uses $MRNA to locate the memerna directory.

Run git submodule init and git submodule update to pull in external dependencies.
Memerna requires a modern C++ compiler that supports C++20.

./build.py

Then run from $PREFIX/memerna. No guarantees this runs or even builds on Windows.

### Running include-what-you-use
./build.py --iwyu --no-build

Then from the cmake build directory:
make -j$(nproc) 2> /tmp/iwyu.out

Then:
iwyu-fix-includes --nocomments --blank_lines --nosafe_headers < /tmp/iwyu.out

### Running the tests
Run from $MRNA/run_tests after building.

### Fuzzing

#### Randomized fuzzing
```
make -j32 && ./fuzz -rd $MRNA/extern/miles_rnastructure/data_tables/ \
  6 8 --print-interval 5 --no-partition --no-subopt
```

Use the -no-table-check option to only compare the result of memerna vs another
program, rather than the internal dp tables.

#### AFL
To run AFL++, first build the afl binary with:

./build.py --type relwithdebinfo --compiler afl-fast

```
sudo sh -c 'echo core >/proc/sys/kernel/core_pattern'
AFL_AUTORESUME=1 AFL_IMPORT_FIRST=1 AFL_TESTCACHE_SIZE=500 AFL_SKIP_CPUFREQ=1 \
  afl-fuzz -x $MRNA/extern/afl/fuzz/dict.dct -m 2000 -t 2000 \
  -i $MRNA/extern/afl/fuzz/testcases -o ./afl -- ./fuzz --afl \
  -rd $MRNA/extern/miles_rnastructure/data_tables/ \
  --mfe-rnastructure --no-subopt --no-partition
```

Reproducing a crash:
```
cat ./afl/default/crashes/<crash>  | ./fuzz --afl ...
```

Minimising test cases:
```
afl-tmin -i case -o ./afl/min -- ./fuzz --afl \
  -rd $MRNA/extern/miles_rnastructure/data_tables/ \
  -r --table-check
```

### Useful commands

./scripts/run_benchmarks.py
Runs benchmarking for various packages.

./scripts/run.py
Supports EFN or folding through -e and -f respectively. Benchmarking through -b.
Can also specify memevault RNA via -mv <memevault name>.

./scripts/parse_data.py
Parses the original data files from orig_data and outputs them in memerna format in data/.

./scripts/extract_subsequence.py
Gives you the subsequence (0 indexed) between st, en of an RNA.

./scripts/crop.py
Crops images in a directory.

./scripts/convert.py
Converts between db and ct files.
