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
  year = {2025},
  month = {7},
  title = {{memerna Software Package}},
  url = {https://github.com/Edgeworth/memerna},
  version = {0.2.0},
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
  publisher={Elsevier},
  doi={10.1016/j.jmb.2024.168819}
}
```

## Building (Ubuntu 22.04 LTS)

Ubuntu 22.04 LTS is the oldest supported Ubuntu release. Its stock Boost 1.74
does not contain Boost.Describe, so the setup below installs Boost 1.83 from
`ppa:mhier/libboost-latest`. On newer distributions, just having a C++ compiler
and cmake is likely sufficient.

The [Dockerfile](Dockerfile) uses vanilla Ubuntu 22.04 and commands that also
work on bare metal. Build and test the image, or run an installed program, with:

```sh
just docker-test
docker run --rm memerna:ubuntu-min fold GCGACCGGGGCUGGCUUGGUAA
```

The image also contains `efn`, `fuzz`, `harness`, `partition`, and `subopt`.

### Environment setup

Install the following packages:

```sh
sudo apt install build-essential cmake git libmpfr-dev software-properties-common
sudo add-apt-repository ppa:mhier/libboost-latest
sudo apt install libboost1.83-dev
```

A direct CMake build does not require Python. The Python tools require uv and
use Python 3.14; uv downloads the interpreter when needed. The repository
command recipes require just. The Dockerfile contains the pinned installation
commands and the additional test dependencies.

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

Alternatively, memerna can be compiled using the rnapy python helper script, if
you have uv installed. This will output builds into `<prefix>/memerna/<config>`,
where the prefix defaults to $HOME/bin.

```sh
uv run python -m rnapy.run build --kind release
```

Optionally, you may set `CPM_SOURCE_CACHE` to a directory to cache external
dependencies. This can be useful if you are building memerna with many separate
configurations.

#### Compiler and toolchain support

GCC 11 or newer is supported. Recent Clang releases should work. Clang 16 might
also work, but Clang builds are not guaranteed. The Ubuntu 22.04 image tests
GCC 11.

### Build configuration

memerna has some build options that can be set with cmake (or passed via the
rnapy helper build script) that may be useful for some users. See the cmake
configuration for a full explanation.

- USE_MPFR: use MPFR for arbitrary precision floating point
- ENABLE_LOGGING: enable project logging
- INDEX_BITS: the number of bits to use for indices and sequence lengths
- ENERGY_PRECISION: the number of decimal places to use for energy calculations
- FLOAT_PRECISION: the number of significant digits to use for floats
- USE_NATIVE_ARCH (cmake only): optimise release builds for the build machine's
  CPU. It is on by default; turn it off for binaries that run on other machines.

ENERGY_PRECISION is by default 2, meaning energies are specified to two decimal
places. If you only need 1 decimal place, you can set this to 1 and performance
will be better.

## Running

By default, memerna auto-selects the backend and algorithm for the requested
energy model and configuration. Use `--backend base`, `--backend canon`, or
`--backend stack` only when you want to force a backend, and use `--dp-alg`,
`--subopt-alg`, or `--pfn-alg` only when you want to force a specific algorithm.

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

With the default energy and trace settings, T04/T12 auto-selection chooses the
sparse base algorithm and T22 auto-selection chooses the stack algorithm.

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

The low-memory variants evict cached expansion vectors and compute them again
after cache misses.

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
./build/fuzz --mfe --mfe-table --subopt --pfn --random-seeds --energy-model t04 \
  1 30
```

Exhaustive fuzzing:

```bash
./build/fuzz --mfe --mfe-table --subopt --pfn --energy-model t04 \
  --enumerate 1 10
```

T22 MFE and suboptimal fuzzing:

```bash
./build/fuzz --mfe --mfe-table --subopt --random-seeds --random-pf \
  --energy-model t22 1 200
```

Random fuzzing continues until it finds a difference or is interrupted. Use
`--fuzz-time-secs N` to stop after N seconds.

`./build/fuzz` loads all supported backends by default. Use `--backends` only to
restrict the comparison set, for example when comparing against an external
package that needs a single memerna backend.

See [CODING.md](CODING.md) for developer fuzzing workflows and implementation
notes.

## License

For any commercial applications of this software, please contact the author for
a license.
