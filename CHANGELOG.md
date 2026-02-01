# Unreleased (since 0.2.1)

## 2026/02
- Refactored configuration passing: Remove CtxCfg and EnergyCfg from Model, pass as function parameters
- Backend can now load dynamically based on algorithm
- Renamed data_dir + seed to variant data_src

## 2026/01
- Added CLAUDE.md with codebase instructions
- Renamed PseudofreeModel to PseudofreeCfg, moved to API
- Standardised backend algorithm calling conventions
- PseudofreeCfg passed as parameter instead of stored in model

---

# 0.2.1 (2025-12-27)

## 2025/12
- Added citation to README
- Fixed bugs and stopped checking in benchmark.json

## 2025/09
- Fixed spdlog corrupting stdout

## 2025/08
- Added explicit pseudofree energy specification to fuzzing

## 2025/07
- Work on multiloop_c support (ongoing)

## 2025/06
- Added versioning support with version.h
- Added multiloop closure cost (multiloop_c) to energy model and MFE
- Refactored Python code structure
- Migrated from pandas to polars for analysis
- Added py.typed file for type checking
- Updated runner output formats

## 2025/05
- Improved subopt perf benchmarking with restart support
- Fixed lowmem mode
- Added CPU affinity support for benchmarking
- Various Python harness improvements

## 2025/04
- Added pseudofree energy support to partition function (pfn) opt algorithm
- Fixed partition function bugs and improved floating point behavior

## 2025/03
- Implemented N2 cache for suboptimal folding (new expansion_cache and lru_vector_cache)
- Implemented pseudofree energy for base partition function (pfn)
- Implemented pseudofree energy for base suboptimal algorithms

## 2025/01
- **Implemented D2 (dangling ends model)** for base backend
- Implemented pseudofree energy for base sparse MFE algorithm
- Added pseudofree energy support to Lyngso sparse algorithm
- Implemented pseudofree energy for hairpin precomputation

---

# 0.2.0 (2025-06-09)

## 2024/09
- Added pseudofree energy support to base energy function (efn)
- Work on adding SHAPE/pseudofree to base model
- Added extra partition tables
- Added subopt partition function fuzzing
- Added partition function epsilon parameter

## 2024/06
- **Large refactor introducing backends architecture**
  - Three backends: base, baseopt, and stack
  - `base`: Full-featured, supports pseudofree energy, all algorithm variants
  - `baseopt` (default): Fastest, strips pseudofree, requires multiloop_c==0
  - `stack`: Alternative energy model with penultimate stacking (T22), no PFN
  - Reorganized code from models/ to backends/
  - Reorganized data model directories (t04p1 -> t04-p1-base, etc.)
  - Tests reorganized to src/tests/

## 2024/05
- Added subopt performance runners
- Added delta calculation for subopt perf runner
- Added count-only mode for subopts
- Moved experiments out of main repo

## 2024/04
- Added SparseRNAFolD to rnapy harness
- Added runner for memerna 0.1 comparison
- Added limits to fold perf runner

---

# 0.1 (2024-01-04)

## 2024/01
- **Implemented persistent data structure for suboptimal folding** (t04 and t22)
  - Uses static cache for 2x performance improvement
  - Optimized priority queue operations
- Improved random seeding for fuzzing
- Updated algorithm names

## 2023/11
- Work on persistent data structures for suboptimal folding
- Updated subopt code structure
- Moved expansion struct
- Python code reformatting and fixes

## 2023/10
- Updated CPM package manager

## 2023/08
- **Dynamic algorithm selection** - basic implementation
- **CPM package manager integration** for dependencies (fmt, spdlog, etc.)
- **Time-based option for suboptimal folding** - stop after time limit
- **No-coax and no-ctd support for t22**
- Added energy config support structure
- Fixed integer overflow for pseudofree tests
- Fixed MPFR for latest fmtlib

## 2023/07
- Fixed `rnapy.run build --regenerate` failing when build folder doesn't exist

---

# 2023/06

Added suboptimal folding for t22:
- Uses same Expansion struct as for randomised traceback
- Pairs are handled explicitly in the State rather than implicitly like in t04
  by looking at the DP tables (for DP_P)

Added SHAPE support for t22:
- add pf_paired and pf_unpaired vectors to the energy model - these are
  pseudofree energy penalties for each nucleotide
- Each substructure handles paired pseudofree energy for its opening pair.
- Closing pairs are the opening pair of the next substructure, and are handled
  by it.
- Unlike RNAstructure:
  - unpaired pseudofree energy is applied for special hairpins
    hairpins.
  - paired pseudofree energy is not double counted inside helices or missed for
    lonely pairs
  - paired pseudofree energy is counted across a size 1 bulge loop
