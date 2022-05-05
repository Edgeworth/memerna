# memerna

Check out docs/thesis/thesis.pdf for my thesis which will explain a bit about
what this is all about.

In all cases where an ordering of base_t p is used (e.g. data tables), it will be ACGU.

### Building

The rest of this document uses $MRNA to locate the memerna directory.

Run git submodule init and git submodule update to pull in external dependencies.
Memerna requires a modern C++ compiler that supports C++20.

python -m rnapy.run build

Then run from $PREFIX/memerna. No guarantees this runs or even builds on Windows.

### Running include-what-you-use

python -m rnapy.run build --iwyu --no-build

Then from the cmake build directory:
make -j$(nproc) 2> /tmp/iwyu.out

Then:
iwyu-fix-includes --nocomments --blank_lines --nosafe_headers < /tmp/iwyu.out

### Running the tests

Run from $MRNA/run_tests after building.

### Fuzzing

```
python -m rnapy.run afl-fuzz --help
```

Reproducing a crash:

```
cat ./afl/default/crashes/<crash>  | ./fuzz --afl ...
```

Minimising test cases:

```
python -m rnapy.run fuzz-min <crash-file>

```
