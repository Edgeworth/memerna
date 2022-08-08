alias t := test
alias f := fix
alias u := update

default:
  @just --list

test:
  # Check various build configurations and run tests. Good to run before push.
  # Uses :::+ to match mpfr-rnastructure and no-mpfr-no-rnastructure
  # Necessary to compile with multiple compilers, as some issues have only shown
  # themselves on a specfic compiler.
  parallel --progress --halt soon,fail=1 --jobs 8 python -m rnapy.run build --test {} ">" /dev/null ::: \
    --kind=debug --kind=relwithdebinfo ::: --sanitizer=asan \
    --sanitizer=tsan --sanitizer=ubsan ::: \
    --rnastructure --no-rnastructure :::+ --mpfr --no-mpfr ::: \
    --compiler=clang --compiler=default

fix:
  SETUPTOOLS_USE_DISTUTILS=stdlib pre-commit run --all-files

update:
  pre-commit autoupdate
  poetry update
