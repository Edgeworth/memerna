# Unreleased (on main, since 0.2.0)

## 2026/02
- Large refactors

## 2025/07
- Added explicit pseudofree energy specification to fuzzing
- Started adding multiloop_c to sparse MFE and subopt algorithms

## 2025/06
- Added versioning support with version.h
- Added multiloop closure cost (multiloop_c) to energy model and MFE
- Migrated from pandas to polars for analysis

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

# Released 0.2.1 (2025-12-27)

- Fixed spdlog corrupting stdout
- Fixed missing newline in seq file output
- Fixed wrong lowmem mode for persistent subopt algorithm

# Released 0.2.0 (2025-06-09)

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
    - unpaired pseudofree energy is applied for special hairpins hairpins.
    - paired pseudofree energy is not double counted inside helices or missed
      for lonely pairs
    - paired pseudofree energy is counted across a size 1 bulge loop
