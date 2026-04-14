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
    --no-logging --logging :::+ \
    --rnastructure --no-rnastructure :::+ --mpfr --no-mpfr ::: \
    --compiler=clang --compiler=default --compiler=afl-fast ::: \
    --index-bits=8 --index-bits=16 --index-bits=32

bench:
  # Run benchmarks.
  poetry run python -m rnapy.run build --bench --bench-output \
    ./benchmark.json --kind=release

fuzz *args:
    #!/usr/bin/env bash
    set -euo pipefail
    if [[ $# -lt 1 ]]; then
      echo "usage: just fuzz <fuzz_exec> [--time-secs <secs>] [--no-pfn]"
      exit 1
    fi
    # `set positional-arguments` passes recipe args to the shebang interpreter.
    # `parallel --shebang` treats that extra argv as an input file, so run
    # `parallel` explicitly from bash instead.
    fuzz_exec=""
    fuzz_time_secs=600
    include_pfn=true
    while [[ $# -gt 0 ]]; do
      case "$1" in
        --time-secs)
          if [[ $# -lt 2 ]]; then
            echo "missing argument for --time-secs"
            exit 1
          fi
          fuzz_time_secs="$2"
          shift
          ;;
        --no-pfn)
          include_pfn=false
          ;;
        --*)
          echo "unknown fuzz arg: $1"
          exit 1
          ;;
        *)
          if [[ -n "$fuzz_exec" ]]; then
            echo "unexpected positional fuzz arg: $1"
            exit 1
          fi
          fuzz_exec="$1"
          ;;
      esac
      shift
    done
    if [[ -z "$fuzz_exec" ]]; then
      echo "missing fuzz_exec"
      exit 1
    fi
    printf -v fuzz_exec_q '%q' "$fuzz_exec"
    base_cfgs=(
      "--mfe --mfe-table 1 100"
      "--mfe --mfe-table --subopt 1 100"
    )
    if [[ "$include_pfn" == true ]]; then
      base_cfgs+=("--mfe --mfe-table --subopt --pfn 1 30")
    fi
    energy_models=(t04 t12 t22)
    ctds=(none d2 no-coax all)
    random_model_ranges=(
      "--random-min-energy -1.0 --random-max-energy 1.0"
      "--random-min-energy -1.0 --random-max-energy 0.0"
      "--random-min-energy 0.0 --random-max-energy 1.0"
      "--random-min-energy 0.0 --random-max-energy 0.0"
    )
    random_pf_ranges=(
      ""
      "--random-pf --random-pf-min-energy -10 --random-pf-max-energy 10"
      "--random-pf --random-pf-min-energy -1 --random-pf-max-energy 1"
      "--random-pf --random-pf-min-energy -10 --random-pf-max-energy -0.1"
      "--random-pf --random-pf-min-energy -1 --random-pf-max-energy -0.1"
      "--random-pf --random-pf-min-energy 0.1 --random-pf-max-energy 10"
      "--random-pf --random-pf-min-energy 0.1 --random-pf-max-energy 1"
    )
    cmds=()
    for base_cfg in "${base_cfgs[@]}"; do
      for energy_model in "${energy_models[@]}"; do
        for ctd in "${ctds[@]}"; do
          for random_model_range in "${random_model_ranges[@]}"; do
            for random_pf_range in "${random_pf_ranges[@]}"; do
              cmd="SPDLOG_LEVEL=err ${fuzz_exec_q} ${base_cfg} --energy-model ${energy_model} --ctd ${ctd} --random-seeds"
              cmd+=" ${random_model_range} --fuzz-time-secs ${fuzz_time_secs}"
              if [[ -n "${random_pf_range}" ]]; then cmd+=" ${random_pf_range}"; fi
              cmds+=("${cmd}")
            done
          done
        done
      done
    done
    echo "Running ${#cmds[@]} fuzzers for ${fuzz_time_secs}s each with $(nproc) jobs"
    parallel --ungroup --halt now,fail=1 --jobs "$(nproc)" bash -lc ::: "${cmds[@]}"

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
