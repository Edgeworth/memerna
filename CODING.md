# Style notes

Set up pre-commit.

## C++

- C++20
- Format according to the .clang-format
- Prefer composition and data oriented APIs. Use inheritance carefully, if at
  all.
- Prefer non-const pointer for output parameters, const-ref for params,
  const-pointer for optional references.

## Directories

- cmake: CMake scripts
- extern: external projects and data
- data: energy model data for memerna
- rnapy: Python code
- rnapy_tests: Python tests
- src: C++ source and tests

### memerna API notes

There are two conceptual APIs:

- A high level API, using Ctx, which makes it easy to do regular RNA
  operations, such as MFE folding with the Turner 2004 model. These
  bundle the outputs into \*Result structs.
- A low level API, where each thing does one thing. These take their inputs as
  basic building block data structures and shouldn't take high level API types.

If a submodule needs some configuration, it should declare a \*Cfg struct, e.g.
BackendCfg, EnergyCfg, SuboptCfg, etc. If it wants to be used from the command
line it should also declare Opt options and a FromArgParse static method on the
config struct.

In general, Cfg structs should not include other Cfg structs. If code needs
multiple Cfg structs, take them separately. This prevents
structs from being included multiple times and potentially conflicting. Maybe
an exception if it's a leaf Cfg struct that won't be included in anything else.

The Cfg file should declare a RegisterOpts function which registers options for
that Cfg struct. It should also call RegisterOpts on any other modules it
would need, if it wanted to be invoked on the command line. e.g. Subopt
should call RegisterOpts for energy and model.

Some notes on the code structure:

## rnapy directory

Contains Python code. Organised as follows:

- analysis: Library for analysing data generated from RNA programs.
- bridge: Library for running various RNA programs.
- build: Library for building and running memerna programs.
- data: Library for getting and manipulating RNA data.
- model: Python RNA and configuration types.
- programs: Programs invoking various scripts functionality.
- util: Misc common utilities.

## Test directories

C++ tests are in `src/tests`, and Python tests are in `rnapy_tests`. Their
structure should mirror the source directories.

## Development checks

`just test` runs the Python tests, then builds and tests representative C++
configurations. These environment variables control its resource use:

- `MRNA_TEST_JOBS`: number of concurrent configurations
- `MRNA_TEST_JOB_MEMORY`: estimated memory per configuration used to select and
  suspend concurrent work, default `4G`
- `MRNA_BUILD_JOBS`: compiler jobs within one configuration
- `MRNA_CLEAN_BUILDS=1`: use a temporary build prefix for each configuration

Each C++ configuration reports elapsed time and maximum resident memory.

To run include-what-you-use, configure a build with:

```sh
uv run python -m rnapy.run build --iwyu --no-build
```

Then run these commands from the generated CMake build directory:

```sh
cmake --build . --parallel "$(nproc)" 2> /tmp/iwyu.out
iwyu-fix-includes --nocomments --blank_lines --nosafe_headers < /tmp/iwyu.out
```

## Fuzzing

`just fuzz` runs the selected fuzz executable across the supported energy-model
and CTD configurations. For example:

```sh
just fuzz ./build/fuzz --time-secs 600
```

Use the Python command for AFL++ setup and execution:

```sh
uv run python -m rnapy.run afl-fuzz --help
```

The Python tools accept these path environment variables:

- `MRNA`: memerna source directory
- `MEMEVAULT`: memevault database file
- `MRNA_DIST`: memerna build directory
- `MRNA01_DIST`: memerna 0.1 build directory
- `RNASTRUCTURE`: RNAstructure source directory with an in-tree build
- `SPARSEMFEFOLD`: SparseMFEFold directory
- `SPARSERNAFOLD`: SparseRNAFolD directory
- `VIENNARNA`: ViennaRNA directory
- `LINEARFOLD`: LinearFold directory

## Energy-model implementation

The `base` and `canon` backends support the T04 and T12 parameter sets. The
`stack` backend also supports T22 and loads `penultimate_stacking.data`.

`IsContinuous` treats a base-pair stack and a size-one bulge as continuous. It
does not treat an internal loop as continuous. The T22 penultimate-stack model
uses this distinction.

Hairpin, bulge, and internal-loop initiation values above 30 use the value at
index 30 plus the extrapolation formula in `ModelBase`.

When bases are used as array indices or data-table dimensions, their order is
A, C, G, U (`A = 0`, `C = 1`, `G = 2`, `U = 3`). Do not change these values.

## Energy data tables

`LoadFromModelPath` loads these common files from each model directory:

- `stacking.data` and `terminal.data`
- `hairpin.data` and `hairpin_initiation.data`
- `bulge_initiation.data`
- `internal_initiation.data`, `internal_1x1.data`, `internal_1x2.data`, and
  `internal_2x2.data`
- `internal_2x3_mismatch.data` and `internal_other_mismatch.data`
- `dangle3.data`, `dangle5.data`, and `misc.data`

Special hairpin and internal-loop table values do not include AU/GU closing-pair
penalties. The energy functions add those penalties separately.

## Releasing

Update CHANGELOG.md and the version numbers in CMakeLists.txt and pyproject.toml,
then run `uv lock`.

Create a new branch from main - branch has major.minor version number.

```bash
git checkout -b release/0.2
```

Add a tag for the release, using the format `release-tag/0.2.0` - it contains
the full version number, including patch.

```bash
git tag release-tag/0.2.0
git push origin release/0.2
git push origin release-tag/0.2.0
```

Then create a release on GitHub, using the tag `release-tag/0.2.0` and the
title "v0.2.0" and description "Release 0.2.0".
