# memerna

In all cases where an ordering of base_t p is used (e.g. data tables), it will be ACGU.


### Useful commands

./scripts/run_efn.py --rnastructure-loc=<> -f <ct or db filename>
Runs EFN for various packages. Can also specify memevault RNA via -k <memevault name>.

./scripts/run_fold.py --rnastructure-loc=<> -f <ct or db filename>
Runs folding for various packages

./scripts/run_benchmarks.py
Runs benchmarking for various packages.

./scripts/run.py
Supports EFN or folding through -e and -p respectively. Benchmarking through -b.

./scripts/plot_benchmarks.py
Plots the output of the benchmarking.

./scripts/parse_data.py
Parses the original data files from orig_data and outputs them in memerna format in data/.

./scripts/extract_subsequence.py
Gives you the subsequence (0 indexed) between st, en of an RNA.

./scripts/crop.py
Crops images in a directory.

./scripts/convert.py
Converts between db and ct files.
