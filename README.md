# memerna

## TODO:

- higher precision energy
  use decimal point format?
  need to update python, programs - use functions to convert
  find all random energy constants, wrap in function to convert.
  search for 100
  allow CAP string
  allow MAX string
- update data tables with higher precision energy, compare to RNAstructure
- ComputeTables0,1,2,3 => better names

## s22 notes (impl)

split code into t04 model and others

compute: code for doing computations with models
compute/t04: general code for t04 model
compute/brute: brute force evaluation glue code
compute/brute/t04: t04 model MFE computation
compute/energy: glue code for energy computations
compute/energy/t04: energy computations in the t04 model
compute/mfe: glue code for MFE computation
compute/mfe/t04: t04 model MFE computation
compute/partition: glue code for partition function computation
compute/partition/t04: t04 model partition function computation
compute/subopt: glue code for suboptimal folding computation
compute/subopt/t04: t04 model suboptimal folding computation
compute/traceback: glue code for traceback computation
compute/traceback/t04: t04 model traceback computation
fuzz: glue code for fuzzing
fuzz/t04: t04 fuzzing

One Cfg for all models. Maybe later need model specific config.

Ctx:

EnergyModel:

same for subopt, partition etc

# s22 notes

2004 model:

- Special stacks of length > 2 base pairs are not handled.
- Internal loops are limited to size 30 in most implementations in memerna.

2012 model: GU penalty removed from AU/GU penalty
2022 model:

- AU penalty removed as well
- stacking params updated
- sequence dependent parameters for terminal base pairs based on the penultimate
  pair

To ask Dave:
need to apply to all loop types - affects energy model, exterior loop, mfe.

- hairpin loop
  - Special hairpins should remain unchanged?
  - Special GU closure (with GG) should remain?
- bulge loop
- internal loop
  - special internal loops should remain unchanged?
  - AU/GU penalty is different for internal loops, but make it the new parameters too?
- multiloop branches
  - hack values for multiloop coefficients should stay the same?
- Still using 10\*kcal/mol. Should switch precision? What accuracy changes did
  Dave see?
- check initiation/symmetry values
  energy model:

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
