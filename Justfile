alias t := test
alias f := fix
alias u := update
set positional-arguments

default:
  @just --list

test:
  poetry run pytest
  # Check various build configurations and run tests. Good to run before push.
  # Uses :::+ to match mpfr-rnastructure and no-mpfr-no-rnastructure etc to
  # save some compilations.
  # Necessary to compile with multiple compilers, as some issues have only shown
  # themselves on a specfic compiler.
  # Use setarch -R to disable ASLR, which can cause issues with thread sanitizer.
  # See https://github.com/google/sanitizers/issues/1716.
  parallel --progress --halt soon,fail=1 --jobs $(nproc) setarch -R poetry run python -m rnapy.run \
    build --test {} ">" /dev/null ::: \
    --kind=debug --kind=relwithdebinfo ::: --sanitizer=asan \
    --sanitizer=tsan --sanitizer=ubsan ::: --float-precision=15 --float-precision=18 :::+ \
    --energy-precision=1 --energy-precision=2 :::+ \
    --rnastructure --no-rnastructure :::+ --mpfr --no-mpfr ::: \
    --compiler=clang --compiler=default --compiler=afl-fast

bench:
  # Run benchmarks.
  poetry run python -m rnapy.run build --bench --bench-output \
    ./benchmark.json --kind=release

fuzz $fuzz_exec:
  #!/usr/bin/env bash
  set -euo pipefail
  # `set positional-arguments` passes recipe args to the shebang interpreter.
  # `parallel --shebang` treats that extra argv as an input file, so run
  # `parallel` explicitly from bash instead.
  printf -v fuzz_exec_q '%q' "$fuzz_exec"
  parallel --ungroup --verbose --halt soon,fail=1 bash -lc ::: \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --pfn --energy-model t04 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --pfn --energy-model t04 1 200" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --random-models 1 200" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --pfn --random-models 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --pfn --random-models --ctd none 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --pfn --random-models --ctd d2 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --pfn --random-models --ctd no-coax 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --pfn --random-models --ctd all 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --pfn --random-models --random-pf --ctd none 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --pfn --random-models --random-pf --ctd d2 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --pfn --random-models --random-pf --ctd no-coax 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --pfn --random-models --random-pf --ctd all 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --random-models --ctd none 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --random-models --ctd d2 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --random-models --ctd no-coax 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --random-models --ctd all 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --random-models --ctd none 1 200" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --random-models --ctd d2 1 200" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --random-models --ctd no-coax 1 200" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --random-models --ctd all 1 200" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --random-models --random-pf --ctd none 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --random-models --random-pf --ctd d2 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --random-models --random-pf --ctd no-coax 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --random-models --random-pf --ctd all 1 30" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --random-models --random-pf --ctd none 1 200" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --random-models --random-pf --ctd d2 1 200" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --random-models --random-pf --ctd no-coax 1 200" \
    "${fuzz_exec_q} --mfe --mfe-table --subopt --random-models --random-pf --ctd all 1 200"

afl-setup:
  #!/usr/bin/env bash
  if [[ "$(cat /proc/sys/kernel/core_pattern)" != "core" ]]; then
    echo core | sudo tee /proc/sys/kernel/core_pattern
  else
    echo "core_pattern already set"
  fi

# Minimize all unique AFL crash files.
afl-tmin afl_dir *args:
  #!/usr/bin/env bash
  afl_dir="$1"
  shift
  mapfile -t crashes < <(find "$afl_dir"/ -path '*/crashes/*' -type f -not -name 'README.txt' \
    -exec md5sum {} + | sort | uniq -w 32 | awk '{print $2}')
  if [ "${#crashes[@]}" -eq 0 ]; then
    echo "No crashes found."
    exit 0
  fi
  echo "Found ${#crashes[@]} unique crashes."
  poetry run python -m rnapy.run afl-fuzz-min --compiler=afl-fast --kind=relwithdebinfo \
    --mfe --mfe-table --subopt "$@" "${crashes[@]}"

afl-fuzz *args: afl-setup
  poetry run python -m rnapy.run afl-fuzz --compiler=afl-fast --kind=relwithdebinfo \
    --mfe --mfe-table --subopt --random-pf "$@"

afl-run *args:
  poetry run python -m rnapy.run afl-fuzz-run --compiler=afl-fast --kind=relwithdebinfo \
    --mfe --mfe-table --subopt --random-pf "$@"

fix:
  pre-commit run --all-files

update:
  poetry run poetry up --latest
  poetry update

check:
  poetry check
  poetry run mypy --install-types --non-interactive rnapy rnapy_tests
  pre-commit run --all-files
