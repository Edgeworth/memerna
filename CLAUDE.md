# CLAUDE.md

RNA secondary structure prediction library with MFE, suboptimal folding, and
partition function algorithms. Optimized for performance with sparse algorithms.

## Build notes
- Test using `just test`
- Check code using `just check`

## Coding notes

- Parameters to functions that are hard to understand, e.g. nullptr should be
  written with the parameter name: `Function(param1, /*important_param=*/
  nullptr);`
- Exhaustive returning switch cases should not use `default:` and should have
  `unreachable()` after them.
- Do not use braces for single-line if/else bodies, unless at least one branch
  is multi-line.

## Critical Invariants

- `Energy` is a 4-byte aligned int32_t. MAX_E (0x0F0F0F0F) and CAP_E (0x07070707)
  are magic values for memset - never change them.
- ENERGY_PRECISION and FLOAT_PRECISION are compile-time only. Data files are
  precision-specific - mismatches will cause incorrect results.
- `baseopt` backend enforces `multiloop_c == ZERO_E`. `base` and `stack`
  support non-zero multiloop closure costs.
- `stack` backend has no partition function support and no sparse algorithms.
- CTDs have dual representation (per-base array and branch deque) that must
  stay synchronized. Per-base stores CTD on right side for outer loops, left
  side for inner multiloops.
- Secondary structure `s[i] = j` implies `s[j] = i` (bidirectional).

## Backend Architecture

Three backends exist with deliberate tradeoffs:
- `base`: Full-featured, supports pseudofree energy, all algorithm variants
- `baseopt` (default): Fastest, strips pseudofree, requires multiloop_c==0
- `stack`: Alternative energy model with penultimate stacking (T22), no PFN

The backends share common code in `backends/common/` but each implements its
own DP state. `stack` uses a 3-array variant system (nostack, penult, standard)
while `base`/`baseopt` use a single 6-type index system.

## Algorithms

MFE: debug (reference), opt (O(N^3)), sparse-opt (O(N^2) expected). Only
`base`/`baseopt` have sparse.

## Energy Models

T04/T12/T22 are Turner model versions. Key differences:
- T22 adds penultimate stacking (only in `stack` backend)
- T12/T22 remove AU/GU penalties progressively
- Lonely pairs handling differs: OFF/HEURISTIC/ON modes

## Testing Strategy

Fuzzing compares multiple backends against each other (consensus verification).
When RNAstructure is compiled in, it serves as ground truth. DP table comparison
mode (`--mfe-table`) catches bugs where final energy is correct by accident.

## Python (rnapy)

Bridge system uses capability matrix - each package declares unsupported ops.
Constraints checked at runtime. Energies use Decimal for exact comparison.

## If something seems wrong in this file, update it.
