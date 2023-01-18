# memerna

## TODO:

when implemented go through all t22 test cases and check structure with -v

- special internal loops need data tables modified to remove penalty; subtract 0.45 if AU/GU
- need to do the same for t12?

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
  -
- See "Nearest neighbor rules for RNA helix folding thermodynamics: improved end effects"

need to apply to all loop types - affects energy model, efn, exterior loop, mfe.

- augu_penalty, internal_augu_penalty goes away;
- add [4][4][4][4] table for terminal base pairs.

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
