# memerna

memerna is a library and set of programs for RNA folding. It includes algorithms
for folding, suboptimal folding, and the partition function. memerna is very
performant, as it uses both fast algorithms and a highly optimized
implementation.

## Building

memerna was tested to build and run in Ubuntu 2022.04 LTS with up to date
packages. The following instructions assume you are using Ubuntu 2022.04 LTS.
However, memerna is mainly developed on Arch Linux, and is likely to work on any
modern Linux distribution provided the toolchain is new enough.

### Environment setup

Install the following packages:

```sh
sudo apt install build-essential cmake git libboost-dev libmpfr-dev
```

On Ubuntu 2022.04 LTS, the following packages are also required since the
toolchain version is too old (note that this is a PPA, so use at your own risk):

```sh
sudo add-apt-repository ppa:ubuntu-toolchain-r/ppa
sudo apt install gcc-12 g++-12
```

Clone the repo:

```sh
git clone https://github.com/Edgeworth/memerna.git
cd memerna
```

Set up git submodules to pull in external dependencies:

```sh
git submodule update --init --recursive
```

### Compiling

memerna can be compiled directly with cmake. The following commands will build
it for Ubuntu 2022.04 LTS:

```sh
CC=gcc-12 CXX=g++-12 cmake -B build -D CMAKE_BUILD_TYPE=Release .
cmake --build build
```

Alternatively, memerna can be compiled using the rnapy python helper script,
if you have Python 3.11+ and poetry installed. This will output builds
into the given prefix directory, which defaults to $HOME/bin/memerna.

```sh
poetry run python -m rnapy.run build --kind release
```

Optionally, you may set `CPM_SOURCE_CACHE` to a directory to cache external
dependencies. This can be useful if you are building memerna with many separate
configurations.

#### Compiler and toolchain support

| Compiler | Version | Supported |
| -------- | ------- | --------- |
| GCC      | <= 11   | ❌        |
| GCC      | 12      | ✅        |
| GCC      | 13      | ✅        |
| Clang    | <= 15   | ❌        |
| Clang    | 16      | ✅        |

Note that clang 14 and 15 will work with a sufficiently modern standard C++
library (but not gcc 11's, or libc++ 14 or 15's).

### Build configuration

memerna has some build options that can be set with cmake (or passed via the
rnapy helper build script) that may be useful for some users. See the cmake
configuration for a full explanation.

- USE_MPFR: use MPFR for arbitrary precision floating point
- ENERGY_PRECISION: the number of decimal places to use for energy calculations
- FLOAT_PRECISION: the number of significant digits to use for floats

ENERGY_PRECISION is by default 2, meaning energies are specified to two decimal
places. If you only need 1 decimal place, you can set this to 1 and performance
will be better.

## Running

### MFE folding

memerna supports minimum free energy (MFE) folding.

```sh
./build/fold GCGACCGGGGCUGGCUUGGUAA
```

There are several algorithms for MFE folding, which can be specified like so:

```sh
./build/fold --dp-alg sparse-opt GCGACCGGGGCUGGCUUGGUAA
```

Some of the algorithms are listed here:

| Algorithm  | Expected time | Expected memory | Example runtime (1000 nt) |
| ---------- | ------------- | --------------- | ------------------------- |
| sparse-opt | O(N^2)        | O(N^2)          | 0.18 seconds              |
| opt        | O(N^3)        | O(N^2)          | 1.58 seconds              |
| debug      | O(N^3)        | O(N^2)          | 3.87 seconds              |


The sparse-opt algorithm is the default, and is the fastest algorithm.

### Suboptimal folding

memerna supports suboptimal folding. For example:

```sh
./build/subopt --ctd-output --subopt-delta 6 GCGACCGGGGCUGGCUUGGUAA
./build/subopt --subopt-delta 6 GCGACCGGGGCUGGCUUGGUAA
./build/subopt --subopt-strucs 7 GCGACCGGGGCUGGCUUGGUAA
./build/subopt --subopt-time-secs 2.5 GCGACCGGGGCUGGCUUGGUAA
```

There are several algorithms for suboptimal folding, which can be specified like so:

```sh
./build/subopt --subopt-alg iterative --ctd-output --subopt-delta 6 GCGACCGGGGCUGGCUUGGUAA
```

Some of the algorithms are listed here (where k is the number of structures produced). The example
runtime is based on a 100 nt sequence and 1000000 structures.

| Algorithm  | Expected time | Expected memory | Example runtime |
| ---------- | ------------- | --------------- | --------------- |
| iterative  | O(N^2 + k)    | O(N^3)          | 12.78 seconds   |
| persistent | O(N^2 + kN)   | O(N^3 + k)      | 2.71 seconds    |

The iterative algorithm will be faster and use less memory for longer sequences. It's theoretically
possible to implement it using O(N^2) memory trading off for worse time performance, but this is not
currently implemented.

### Partition function

memerna supports computing the partition function. For example:

```sh
./build/partition GCGACCGGGGCUGGCUUGGUAA
```

There are several algorithms for the partition function, which can be specified like so:

```sh
./build/partition --part-alg opt GCGACCGGGGCUGGCUUGGUAA
```

Some of the algorithms are listed here:

| Algorithm | Expected time | Expected memory | Example runtime (1000 nt) |
| --------- | ------------- | --------------- | ------------------------- |
| opt       | O(N^3)        | O(N^2)          | 9.23 seconds              |
| debug     | O(N^3)        | O(N^2)          | 34.43 seconds             |

### Running the tests

```sh
./build/run_tests
```

### Running include-what-you-use

poetry run python -m rnapy.run build --iwyu --no-build

Then from the cmake build directory:
make -j$(nproc) 2> /tmp/iwyu.out

Then:
iwyu-fix-includes --nocomments --blank_lines --nosafe_headers < /tmp/iwyu.out

### Fuzzing

Fuzzing against RNAstructure

```bash
poetry run python -m rnapy.run build --kind relwithdebinfo --rnastructure --energy-precision 1
# Just MFE:
./build/fuzz -rd $MRNA/extern/rnastructure_bridge/data_tables/ --mfe --mfe-rnastructure --mfe-table 1 200
```

### Fuzzing with afl-fuzz

Note that afl-fast seems to cause broken behaviour recently, compared to afl-lto.

```bash
poetry run python -m rnapy.run afl-fuzz --help
```

For example, try this command line:

```bash
poetry run python -m rnapy.run afl-fuzz --kind release --compiler afl-lto \
 --mfe --num-procs 1 --max-len 500 --energy-model t22p2 --seed 1234
```

To fuzz everything:

```bash
poetry run python -m rnapy.run afl-fuzz --kind release --compiler afl-lto \
 --mfe --mfe-rnastructure --mfe-table --part --part-rnastructure \
 --subopt --subopt-rnastructure --num-procs 28 --max-len 100 \
 --energy-model t04p2 --energy-precision 1
```

Fuzz memerna only, with faster settings:

```bash
poetry run python -m rnapy.run afl-fuzz --kind release --compiler afl-lto \
 --mfe --mfe-table --subopt --num-procs 28 --max-len 100 --seed 1234 \
 --energy-model t04p2 --subopt-strucs 100 --subopt-delta 0.2 \
 --energy-precision 2
```

Checking progress:

```bash
afl-whatsup -s $PREFIX/memerna-afl/\*/afl
```

Reproducing a crash:

```bash
cat ./afl/default/crashes/<crash> | ./fuzz --afl ...
```

Minimising test cases:

```bash
poetry run python -m rnapy.run afl-fuzz-min <crash-file>
```

### Environment variables for rnapy

Optionally, it can be useful to set the following variables, for example in
a .env file (used by rnapy):

```bash
MRNA=<set to memerna source directory>
MEMEVAULT=${MRNA}/rnapy/data/memevault.db
MRNA_DIST=${HOME}/bin/memerna/relwithdebinfo-default-64-rnastructure
RNASTRUCTURE=${HOME}/...
SPARSEMFEFOLD=${HOME}/...
VIENNARNA=${HOME}/...
LINEARFOLD=${HOME}/...
```

## Model notes

Turner 1999 model (not implemented)

Turner 2004 model (t04p1, t04p2):

- Adds coaxial stacking
- Lonely pairs are "soft disallowed" - only consider pairs where at least one of
  the adjacent two pairs could be made.
- Updates parameters for terminal mismatches, hairpin, bulge, internal, and multiloops.
- Special stacks of length > 2 base pairs are not handled.
- Internal loops are limited to size 30 in most implementations in memerna.

Update in 2012 model (t12p2):

- GU penalty removed from AU/GU penalty
- GU penalty removed from special hairpins, if closed by GU (none like this)
- Removed 0.45 portion of GU penalty from internal loop GU penalty.
- Stacking parameters changed (for GU stacks, not WC)
- See "Testing the nearest neighbor model for Canonical RNA base pairs" paper

Update in 2022 model (t22p2):

- AU penalty removed as well
- AU penalty removed from special hairpins, if closed by AU
- Lonely pairs implemented fully correctly.
- Removed 0.5 (not 0.45) portion of AU penalty from internal loop AU penalty.
- Stacking parameters changed
- Sequence dependent parameters for terminal base pairs based on the penultimate
  pair (penultimate_stacking.data)
  - Bulge loops of size 1 are continuous and don't have penultimate stacking
    applied inside (same as AU/GU penalty rules).
  - Coaxial stacks are not treated as continuous for penultimate stacking (same
    as AU/GU penalty rules).
- See "Nearest neighbor rules for RNA helix folding thermodynamics: improved end effects"

## Data table notes

- bulge_initiation.data
- dangle3.data
- dangle5.data
- hairpin.data
  Special hairpins. AU/GU penalties NOT baked in.
- hairpin_initiation.data
- internal_1x1.data
  1x1 special internal loops. AU/GU penalties NOT baked in.
  N.B. internal loops are never treated as continuous.
- internal_1x2.data
  2x1 special internal loops. AU/GU penalties NOT baked in.
- internal_2x2.data
  2x2 special internal loops. AU/GU penalties NOT baked in.
- internal_2x3_mismatch.data
  2x3 internal loop terminal mismatch parameters. AU/GU penalties NOT baked in.
- internal_initiation.data
- internal_other_mismatch.data
  Non 2x3 internal loop terminal mismatch parameters. AU/GU penalties NOT baked in.
- misc.data
- stacking.data
- terminal.data

## Misc notes

In all cases where an ordering of base_t p is used (e.g. data tables), it will be ACGU.

## License notes

For any commercial applications of this software, please contact the author for
a license.
