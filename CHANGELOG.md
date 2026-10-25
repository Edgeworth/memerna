# Unreleased (on main, since 0.2.0)

## 2026/08
- Added a docker Ubuntu 22.04 build, test, and runtime image.
- Updated the Python tools to Python 3.14 and replaced mypy with ty.
- Misc fixes.

## 2026/07
- Required CMake 3.16 and Boost 1.77 or newer.
- Misc fixes.

## 2026/06
- Added the canon backend for T04 and T12.
- Misc fixes + optimizations.
- Replaced Poetry with uv for Python dependency management.

## 2026/04
- Removed the baseopt backend and its model data.
- Added compile-time instantiated algorithms (speedup).
- Expanded fuzz coverage and random energy configuration.

## 2026/03
- Fixed sparse candidate pruning when the multiloop closure cost is nonzero.

## 2026/02
- Added automatic backend and algorithm selection.
- Added multiloop closure costs to stack MFE and suboptimal folding.
- Removed the obsolete brute and splay-explorer programs.

## 2025/07
- Added explicit pseudofree energy specification to fuzzing
- Started adding multiloop_c to sparse MFE and subopt algorithms

## 2025/06
- Added multiloop closure cost (multiloop_c) to energy model and MFE
- Migrated from pandas to polars for analysis

# Released 0.2.4 (2026-10-25)

- Turned off sanitizers and RNAstructure in `just test` when their git
  submodules are missing.

# Released 0.2.3 (2026-10-18)

- Added a docker Ubuntu 22.04 build, test, and runtime image.
- Supported GCC 11 and Python 3.10.
- Required Boost 1.77 or newer.
- Disabled RNAstructure by default. Set USE_RNASTRUCTURE=ON to build with it.
- Only required the git submodules for RNAstructure and sanitizer builds.

# Released 0.2.2 (2026-02-13)

- Supported more than 2^31 - 1 suboptimal structures.
- Improved suboptimal performance-runner handling of timeouts and failed
  commands.

# Released 0.2.1 (2025-12-27)

- Fixed spdlog corrupting stdout
- Fixed missing newline in seq file output
- Fixed wrong lowmem mode for persistent subopt algorithm

# Released 0.2.0 (2025-06-09)

## 2025/06
- Added versioning support with version.h

## 2025/05
- Fixed lowmem mode for suboptimal folding
- Added CPU affinity support for benchmarking

## 2025/04
- Added pseudofree energy support to partition function opt algorithm
- Fixed partition function unpaired nucleotide handling

## 2025/03
- Implemented N^2 cache for suboptimal folding (new expansion_cache and lru_vector_cache)
- Implemented pseudofree energy for base partition function
- Implemented pseudofree energy for base suboptimal algorithms

## 2025/01
- Implemented D2 (dangling ends model) for base backend
- Implemented pseudofree energy for base sparse MFE algorithm
- Added pseudofree energy support to Lyngso sparse algorithm
- Implemented pseudofree energy for hairpin precomputation

## 2024/09
- Added pseudofree energy support to base energy function (efn)

## 2024/06
- Large refactor introducing backends architecture
  - Three backends: base, baseopt, and stack
  - `base`: Full-featured, supports pseudofree energy, all algorithm variants
  - `baseopt` (default): Fastest, strips pseudofree, requires multiloop_c==0
  - `stack`: Alternative energy model with penultimate stacking (T22), no PFN
  - Reorganized code from models/ to backends/
  - Reorganized data model directories (t04p1 -> t04-p1-base, etc.)
  - Tests reorganized to src/tests/

## 2024/05
- Moved experiments out of main repo

# Released 0.1 (2024-01-04)

## 2024/01
- Implemented persistent data structure for suboptimal folding (t04 and t22)

## 2023/11
- Started implementing persistent data structures for suboptimal folding
- Refactored subopt code structure and moved expansion struct

## 2023/08
- Dynamic algorithm selection - basic implementation
- CPM package manager integration for dependencies
- Time-based option for suboptimal folding
- No-coax and no-ctd support for t22

## 2023/07
- Fixed `rnapy.run build --regenerate` failing when build folder doesn't exist

# 2023/06

- Added suboptimal folding for t22:
  - Uses same Expansion struct as for randomised traceback
  - Pairs are handled explicitly in the State rather than implicitly like in t04
    by looking at the DP tables (for DP_P)
- Added SHAPE support for t22:
  - add pf_paired and pf_unpaired vectors to the energy model - these are
    pseudofree energy penalties for each nucleotide
  - Each substructure handles paired pseudofree energy for its opening pair.
  - Closing pairs are the opening pair of the next substructure, and are handled
    by it.
  - Unlike RNAstructure:
    - unpaired pseudofree energy is applied for special hairpins.
    - paired pseudofree energy is not double counted inside helices or missed
      for lonely pairs
    - paired pseudofree energy is counted across a size 1 bulge loop
