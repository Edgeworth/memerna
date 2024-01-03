# memerna

## TODO

- use CPM to handle libraries
- implement d2
  - add to mfe_fastest, see how much it affects perf
  - add all variants to all mfe, traceback, suboptimal, partition function
  - add to energy model
  - add to fuzz
  - need to add to ComputeOptimalCtds
  - ctrl f for flush coax etc

need to implement lonely pairs disabling properly for t22. (check this)

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

## RNAstructure data tables notes

- rna.coaxial
- rna.coaxstack
- rna.cov
- rna.dangle
- rna.dynalignmiscloop
- rna.helix_ends
  Table for terminal end penalties (e.g. AU/GU penalties).
- rna.hexaloop
  Special hairpin loops of length 6. AU/GU penalties baked in.
- rna.int11
  1x1 special internal loops. AU/GU penalties baked in.
- rna.int21
  2x1 special internal loops. AU/GU penalties baked in.
- rna.int22
  2x2 special internal loops. AU/GU penalties baked in.
- rna.loop
- rna.miscloop
- rna.param_map
- rna.stack
  Stacking parameters for helices.
- rna.tloop
  Special hairpin loops of length 4. AU/GU penalties baked in.
- rna.triloop
  Special hairpin loops of length 3. AU/GU penalties baked in.
- rna.tstack
- rna.tstackcoax
- rna.tstackh
  Hairpin loop terminal mismatch parameters.
  AU/GU penalties and other special hairpin parameters baked in.
- rna.tstacki
  Internal loop terminal mismatch parameters. Internal loop AU/GU penalties baked in.
- rna.tstacki1n
  Bulge loop terminal mismatch parameters. AU/GU penalties baked in.
- rna.tstacki23
  2x3 internal loop terminal mismatch parameters. Internal loop AU/GU penalties baked in.
- rna.tstackm

## Misc notes

In all cases where an ordering of base_t p is used (e.g. data tables), it will be ACGU.

### Environment variables

Optionally, it can be useful to set the following variables, for example in
a .env file (used by rnapy):

```bash
MRNA=...
MEMEVAULT=${MRNA}/rnapy/data/memevault.db
MRNA_DIST=${HOME}/bin/memerna/relwithdebinfo-default-64-rnastructure
RNASTRUCTURE=${HOME}/...
SPARSEMFEFOLD=${HOME}/...
VIENNARNA=${HOME}/...
LINEARFOLD=${HOME}/...
```

### Building

The rest of this document uses $MRNA to locate the memerna directory.
Memerna requires a modern C++ compiler that supports C++20.

Set up git submodules to pull in external dependencies.:

```sh
git submodule init
git submodule update
```

Install dependencies: boost, python poetry. Then we can build:

```sh
poetry run python -m rnapy.run build
```

Then run from $PREFIX/memerna. No guarantees this runs or even builds on Windows.

### Compiler and toolchain support

| Compiler | Version      | Supported |
|----------|--------------|-----------|
| GCC      | <= 11        | ❌        |
| GCC      | 12           | ✅        |
| GCC      | 13           | ✅        |
| Clang    | <= 13        | ❓        |
| Clang    | 14           | ✅        |
| Clang    | 15           | ✅        |
| Clang    | 16           | ✅        |

### Running include-what-you-use

poetry run python -m rnapy.run build --iwyu --no-build

Then from the cmake build directory:
make -j$(nproc) 2> /tmp/iwyu.out

Then:
iwyu-fix-includes --nocomments --blank_lines --nosafe_headers < /tmp/iwyu.out

### Running the tests

Run from $MRNA/run_tests after building.

### Fuzzing

Fuzzing against RNAstructure

```bash
poetry run python -m rnapy.run build --kind relwithdebinfo --rnastructure --energy-precision 1
# Just MFE:
./fuzz -rd $MRNA/extern/rnastructure_bridge/data_tables/ --mfe --mfe-rnastructure --mfe-table 1 200
# Partition and subopt are supported, but fuzzing shows differences instantly.
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

## License notes

For any commercial applications of this software, please contact the author for
a license.
