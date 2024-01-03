# memerna

memerna is a library and set of programs for RNA folding. It includes algorithms
for folding, suboptimal folding, and the partition function. memerna is very
performant, as it uses both fast algorithms and a highly optimized
implementation.

This is version 0.1 of memerna. The main branch is actively developed and has
more features (although it is not API or command line compatible with 0.1),
although is slightly less performant. The main branch is recommended in general.

## Building

memerna 0.1 was tested to build and run in Ubuntu 2022.04 LTS with up to date
packages. The following instructions assume you are using Ubuntu 2022.04 LTS.

### Environment setup

Install the following packages:

```sh
sudo apt install build-essential cmake git
```

Clone the repo and checkout v0.1:

```sh
git clone https://github.com/Edgeworth/memerna.git
cd memerna
git switch release/0.1
```

Set up git submodules to pull in external dependencies:

```sh
git submodule update --init --recursive
```

The default configuration builds without MPFR or RNAstructure integration. Most
people will not need either of these. If you want to use MPFR, install the following
dependencies:

```sh
sudo apt install libboost-dev libmpfr-dev
```

### Compiling

memerna can be compiled directly with cmake. The following commands will build
it:

```sh
cmake -B build -D CMAKE_BUILD_TYPE=Release .
cmake --build build
```

There is also the build.py script which wraps cmake for convenience. The
build.py script creates some additional directory structure for each build type.

```sh
./build.py -t release -p ./build
```

#### Compiler and toolchain support

| Compiler | Version      | Supported |
|----------|--------------|-----------|
| GCC      | <= 10        | ❓        |
| GCC      | >= 11        | ✅        |
| Clang    | <= 13        | ❓        |
| Clang    | >= 14        | ✅        |

## Running

The following commands assume memerna is built in the build/ directory and the
current directory is the memerna source directory.

You may pass `--help` to any of the programs to see a list of options.

### MFE folding examples

```sh
./build/fold -memerna-data ./data/ GCGACCGGGGCUGGCUUGGUAA
```

### Suboptimal folding examples

```sh
./build/subopt -memerna-data ./data/ -ctd-output -delta 6 GCGACCGGGGCUGGCUUGGUAA
./build/subopt -memerna-data ./data/ -delta 6 GCGACCGGGGCUGGCUUGGUAA
./build/subopt -memerna-data ./data/ -num 7 GCGACCGGGGCUGGCUUGGUAA
```

### Partition function examples

```sh
./build/partition -memerna-data ./data/ GCGACCGGGGCUGGCUUGGUAA
```

### Running the tests

```sh
./build/run_tests -memerna-data ./data/
```

## Source

- cmake: CMake scripts
- data: energy model data for memerna
- examples: various dot-bracket example folded RNAs
- extern: external projects and data
- scripts: scripts for various things
- src: source
- tests: tests

## License notes

For any commercial applications of this software, please contact the author for
a license.
