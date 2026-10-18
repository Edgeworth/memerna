# memerna

memerna is a library and set of programs for RNA folding. It includes algorithms
for folding, suboptimal folding, and the partition function. memerna is very
performant, as it uses both fast algorithms and a highly optimized
implementation.

## Citation

Please cite the following for the memerna software package:

```
@software{memerna2025,
  author = {Courtney, Eliot and Ward, Max},
  year = {2026},
  month = {10},
  title = {{memerna Software Package}},
  url = {https://github.com/Edgeworth/memerna},
  version = {0.2.3},
}
```

Please cite the following for sparse folding including coaxial stacking:

```
@article{courtney2025memerna,
  title={memerna: Sparse RNA folding including coaxial stacking},
  author={Courtney, Eliot and Datta, Amitava and Mathews, David H and Ward, Max},
  journal={Journal of Molecular Biology},
  volume={437},
  number={3},
  pages={168819},
  year={2025},
  publisher={Elsevier}
}
```

## Building (Ubuntu 22.04 LTS)

Ubuntu 22.04 LTS is the oldest supported Ubuntu release. Its stock Boost 1.74
does not contain Boost.Describe, so the setup below installs Boost 1.83 from
`ppa:mhier/libboost-latest`. On newer distributions, a C++ compiler and cmake
are likely to be sufficient.

The [Dockerfile](Dockerfile) uses vanilla Ubuntu 22.04 and the same tools that
can be installed on bare metal. Build and test the image, or run an installed
program, with:

```sh
just docker-test
docker run --rm memerna:ubuntu-min fold GCGACCGGGGCUGGCUUGGUAA
```

The image also contains `brute`, `efn`, `fuzz`, `harness`, `partition`, and
`subopt`.

### Environment setup

Install the following packages:

```sh
sudo apt install build-essential cmake git libmpfr-dev software-properties-common
sudo add-apt-repository ppa:mhier/libboost-latest
sudo apt install libboost1.83-dev
```

A direct CMake build does not require Python. The Python helper requires Python
3.10 or newer and Poetry. The repository command recipes also require just. The
Dockerfile uses Ubuntu's Python 3.10 and installs Poetry with pip in a virtual
environment.

Clone the repo:

```sh
git clone https://github.com/Edgeworth/memerna.git
cd memerna
```

To build with RNAstructure (USE_RNASTRUCTURE) or sanitizers, set up the git
submodules:

```sh
git submodule update --init --recursive
```

### Compiling

memerna can be compiled directly with cmake. The following commands will build
it for Ubuntu 22.04 LTS:

```sh
cmake -B build -D CMAKE_BUILD_TYPE=Release .
cmake --build build
```

Alternatively, memerna can be compiled using the rnapy python helper script,
if you have Python 3.10+ and poetry installed. This will output builds
into the given prefix directory, which defaults to $HOME/bin/memerna.

```sh
poetry install
poetry run python -m rnapy.run build --kind release
```

Optionally, you may set `CPM_SOURCE_CACHE` to a directory to cache external
dependencies. This can be useful if you are building memerna with many separate
configurations.

#### Compiler and toolchain support

GCC 11 or newer is supported. The Ubuntu 22.04 image tests its default GCC 11
toolchain. Recent Clang releases should work, but Clang builds are not
guaranteed.

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

### Subopt folding

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

Some of the algorithms are listed here (where K is the number of structures produced). The example
runtime is based on a 100 nt sequence and 1000000 structures.

| Algorithm         | Expected time           | Expected memory | Example runtime |
| ----------------- | ----------------------- | --------------- | --------------- |
| iterative         | O(N^3 + KN)             | O(N^3)          | 12.78 seconds   |
| persistent        | O(N^3 + KN + K log K)   | O(N^3 + K)      | 2.71 seconds    |
| iterative-lowmem  | O(N^3 + KN^2)           | O(N^2)          |                 |
| persistent-lowmem | O(N^3 + KN^2 + K log K) | O(N^2 + K)      |                 |

The iterative algorithm will be faster and use less memory for longer sequences. The lowmem variants
use an LRU cache to reduce memory from O(N^3) to O(N^2), at the cost of increased time per structure
from O(N) to O(N^2) due to cache misses.

### Partition function

memerna supports computing the partition function. For example:

```sh
./build/partition GCGACCGGGGCUGGCUUGGUAA
```

There are several algorithms for the partition function, which can be specified like so:

```sh
./build/partition --pfn-alg opt GCGACCGGGGCUGGCUUGGUAA
```

Some of the algorithms are listed here:

| Algorithm | Expected time | Expected memory | Example runtime (1000 nt) |
| --------- | ------------- | --------------- | ------------------------- |
| opt       | O(N^3)        | O(N^2)          | 9.23 seconds              |
| debug     | O(N^3)        | O(N^2)          | 34.43 seconds             |

The partition function is not supported for T22.

### Running the tests

```sh
./build/run_tests
```

## Fuzzing

### Basic fuzzing

General fuzzing:

```bash
./build/fuzz --mfe --mfe-table --subopt --pfn --random-models --energy-model t04 \
  --backends base,baseopt 1 30
```

Exhaustive fuzzing:

```bash
./build/fuzz --mfe --mfe-table --subopt --pfn --energy-model t04 \
  --backends base,baseopt --enumerate 1 10
```

Fuzzing for t22:

```bash
./build/fuzz --mfe --mfe-table --subopt --random-models --random-pf \
  --energy-model t22 --backends stack 1 200
```

### Fuzzing in parallel

The below command runs fuzzing in parallel using GNU parallel with 16 jobs.

```bash
seq 16 | parallel -j 16 -n0 -u './build/fuzz --mfe --mfe-table --subopt --pfn \
  --random-models --energy-model t04 --backends base,baseopt 1 30'
```

### Fuzzing against RNAstructure

```bash
poetry run python -m rnapy.run build --kind relwithdebinfo --rnastructure --energy-precision 1
# Just MFE:
./build/fuzz -rd $MRNA/extern/rnastructure_bridge/data_tables/ --mfe \
  --mfe-rnastructure --mfe-table --energy-model t04 --backends base,baseopt 1 200
```

### Fuzzing with afl-fuzz

Note that afl-fast seems to cause broken behaviour recently, compared to
afl-lto. It's better to use a seeded energy model than random energy models for
afl-fuzz since it will be more stable and easier to debug if it finds something.

```bash
poetry run python -m rnapy.run afl-fuzz --help
```

For example, try this command line:

```bash
poetry run python -m rnapy.run afl-fuzz --kind release --compiler afl-lto \
 --num-procs 16 --mfe --mfe-table --seed 123 --max-len 200 \
 --energy-model t22 --backends stack
```

To fuzz everything:

```bash
poetry run python -m rnapy.run afl-fuzz --kind release --compiler afl-lto \
  --num-procs 16 --energy-precision 1 --mfe --mfe-rnastructure --mfe-table \
  --subopt --max-len 200 --backends base,baseopt
```

Fuzz memerna only, with faster settings:

```bash
poetry run python -m rnapy.run afl-fuzz --kind release --compiler afl-lto \
 --num-procs 16 --mfe --mfe-table --subopt --subopt-strucs 100 \
 --subopt-delta 0.2 --seed 123 --max-len 200 --backends base,baseopt
```

Checking progress:

```bash
afl-whatsup -s $PREFIX/memerna-afl/*/afl
```

Reproducing a crash:

```bash
cat ./afl/default/crashes/<crash> | ./build/fuzz --afl ...
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
MRNA01_DIST=${HOME}/...
RNASTRUCTURE=${HOME}/...
SPARSEMFEFOLD=${HOME}/...
SPARSERNAFOLD=${HOME}/...
VIENNARNA=${HOME}/...
LINEARFOLD=${HOME}/...
```

## Model notes

Turner 1999 model (not implemented)

Turner 2004 model (t04):

- Adds coaxial stacking
- Lonely pairs are "soft disallowed" - only consider pairs where at least one of
  the adjacent two pairs could be made.
- Updates parameters for terminal mismatches, hairpin, bulge, internal, and multiloops.
- Special stacks of length > 2 base pairs are not handled.
- Internal loops are limited to size 30 in most implementations in memerna.

Update in 2012 model (t12):

- GU penalty removed from AU/GU penalty
- GU penalty removed from special hairpins, if closed by GU (none like this)
- Removed 0.45 portion of GU penalty from internal loop GU penalty.
- Stacking parameters changed (for GU stacks, not WC)
- See "Testing the nearest neighbor model for Canonical RNA base pairs" paper

Update in 2022 model (t22):

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
