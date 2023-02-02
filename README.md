# memerna

## TODO:

when implemented go through all t22 test cases and check structure with -v

- special internal loops need data tables modified to remove penalty; subtract 0.45 if AU/GU or 0.5?
- need to do the same for t12?

internal loop / etc au/gu penalties: subtract 0.5 or 0.45?

questions about RNAstructure data tables:
rna.helix_ends.dg: looks like it contains AU/GU penalties. Why?
rna.stack.dg: Why not updated with new stacking parameters?

need to implement lonely pairs disabling properly for t22.

add terminal stacks to:
stacking:
hairpins: what about ((...)) - need terminal stacks here?
internal loops: ((..((...))..))
bulge loops: ((((...)).))

check float bits, partition fn

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
- Removed 0.45 portion of AU penalty from internal loop AU penalty.
- Stacking parameters changed
- Sequence dependent parameters for terminal base pairs based on the penultimate
  pair (penultimate_stacking.data)
  - Bulge loops of size 1 are continuous and don't have penultimate stacking
    applied inside (same as AU/GU penalty rules).
  - Coaxial stacks are not treated as continuous for penultimate stacking (same
    as AU/GU penalty rules).
- See "Nearest neighbor rules for RNA helix folding thermodynamics: improved end effects"

## Data table notes:

- bulge_initiation.data
- dangle3.data
- dangle5.data
- hairpin.data
  Special hairpins. AU/GU penalties baked in.
- hairpin_initiation.data
- internal_1x1.data
  1x1 special internal loops. AU/GU penalties baked in.
  N.B. internal loops are never treated as continuous.
- internal_1x2.data
  2x1 special internal loops. AU/GU penalties baked in.
- internal_2x2.data
  2x2 special internal loops. AU/GU penalties baked in.
- internal_2x3_mismatch.data
  2x3 internal loop terminal mismatch parameters. AU/GU not baked in.
- internal_initiation.data
- internal_other_mismatch.data
  Non 2x3 internal loop terminal mismatch parameters. AU/GU not baked in.
- misc.data
- stacking.data
- terminal.data

## RNAstructure data tables notes

- rna.coaxial
- rna.coaxstack
- rna.cov
- rna.dangle
- rna.dynalignmiscloop
- rna.hexaloop
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
- rna.tloop
- rna.triloop
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

```
MRNA=...
MEMEVAULT=${MRNA}/rnapy/data/memevault.db
MRNA_DIST=${HOME}/bin/memerna/relwithdebinfo-default-64-rnastructure
RNASTRUCTURE=${HOME}/...
SPARSEMFEFOLD=${HOME}/...
VIENNARNA=${HOME}/...
```

### Building

The rest of this document uses $MRNA to locate the memerna directory.

Run git submodule init and git submodule update to pull in external dependencies.
Memerna requires a modern C++ compiler that supports C++20.

Install dependencies: gtest, google/benchmark, fmt

python -m rnapy.run build

Then run from $PREFIX/memerna. No guarantees this runs or even builds on Windows.

### Running include-what-you-use

python -m rnapy.run build --iwyu --no-build

Then from the cmake build directory:
make -j$(nproc) 2> /tmp/iwyu.out

Then:
iwyu-fix-includes --nocomments --blank_lines --nosafe_headers < /tmp/iwyu.out

### Running the tests

Run from $MRNA/run_tests after building.

### Fuzzing

Note that afl-fast seems to cause broken behaviour recently, compared to afl-lto.

```
python -m rnapy.run afl-fuzz --help
```

Reproducing a crash:

```
cat ./afl/default/crashes/<crash>  | ./fuzz --afl ...
```

Minimising test cases:

```
python -m rnapy.run afl-fuzz-min <crash-file>

```
