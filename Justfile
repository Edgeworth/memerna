set positional-arguments
alias t := test
alias f := fix
alias u := update

default:
  @just --list

test:
  # Check various build configurations and run tests. Good to run before push.
  # Uses :::+ to match mpfr-rnastructure and no-mpfr-no-rnastructure etc to
  # save some compilations.
  # Necessary to compile with multiple compilers, as some issues have only shown
  # themselves on a specfic compiler.
  # Use setarch -R to disable ASLR, which can cause issues with thread sanitizer.
  # See https://github.com/google/sanitizers/issues/1716.
  parallel --progress --halt soon,fail=1 --jobs 8 setarch -R poetry run python -m rnapy.run \
    build --test {} ">" /dev/null ::: \
    --kind=debug --kind=relwithdebinfo ::: --sanitizer=asan \
    --sanitizer=tsan --sanitizer=ubsan ::: --float-precision=15 --float-precision=18 :::+ \
    --energy-precision=1 --energy-precision=2 :::+ \
    --rnastructure --no-rnastructure :::+ --mpfr --no-mpfr ::: \
    --compiler=clang --compiler=default

bench:
  # Run benchmarks.
  poetry run python -m rnapy.run build --bench --bench-output \
    ./data/benchmark.json --kind=release

fix:
  pre-commit run --all-files

update:
  poetry run poetry up --latest
  poetry update
  pre-commit autoupdate
