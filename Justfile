alias t := test
alias f := fix
alias u := update
set positional-arguments

default:
  @just --list

test:
  #!/usr/bin/env bash
  set -euo pipefail
  uv run pytest
  test_job_memory="${MRNA_TEST_JOB_MEMORY:-10G}"
  test_jobs="${MRNA_TEST_JOBS:-}"
  read -r -a test_compilers <<< \
    "${MRNA_TEST_COMPILERS:---compiler=clang --compiler=default}"
  if [[ -z "$test_jobs" ]]; then
    # MemAvailable includes reclaimable cache. Leave 25% for other work.
    mem_available_kib="$(awk '$1 == "MemAvailable:" { print $2; exit }' /proc/meminfo)"
    test_job_bytes="$(numfmt --from=iec "$test_job_memory")"
    memory_jobs=$((mem_available_kib * 3 * 1024 / 4 / test_job_bytes))
    if ((memory_jobs < 1)); then memory_jobs=1; fi
    test_jobs="$(nproc)"
    if ((memory_jobs < test_jobs)); then test_jobs="$memory_jobs"; fi
  fi
  echo "Testing build configurations with $test_jobs jobs and $test_job_memory per job"
  # Check representative build configurations. Good to run before push.
  # Uses :::+ to match mpfr-rnastructure and no-mpfr-no-rnastructure etc to
  # save some compilations.
  # Necessary to compile with multiple compilers, as some issues have only shown
  # themselves on a specific compiler.
  # memerna is single threaded, so a thread sanitizer build finds nothing.
  cfgs="$(parallel --keep-order echo ::: \
    --kind=debug --kind=relwithdebinfo ::: --sanitizer=asan \
    --sanitizer=ubsan ::: --float-precision=15 --float-precision=18 :::+ \
    --energy-precision=1 --energy-precision=2 :::+ \
    --no-logging --logging :::+ \
    --rnastructure --no-rnastructure :::+ --mpfr --no-mpfr ::: \
    "${test_compilers[@]}" ::: --index-bits=16 --index-bits=32)"
  base="--kind=debug --float-precision=15 --energy-precision=1 --no-logging"
  base+=" --rnastructure --mpfr"
  for test_compiler in "${test_compilers[@]}"; do
    # Also check one build without a sanitizer, since sanitizers change codegen.
    cfgs+=$'\n'"$base --sanitizer=none --index-bits=16 $test_compiler"
    # 8 bit indices cap the sequence length, so the tests skip their longest
    # sequence. One build covers the size limit without shrinking the matrix.
    cfgs+=$'\n'"$base --sanitizer=ubsan --index-bits=8 $test_compiler"
  done
  # The submodule directories are empty in a source archive and in a clone
  # before git submodule update, so turn off the options that need them.
  if [[ ! -f extern/sanitizers-cmake/cmake/FindSanitizers.cmake ]]; then
    echo "Testing without sanitizers because extern/sanitizers-cmake is missing"
    cfgs="$(sed -E 's/--sanitizer=[a-z]+/--sanitizer=none/' <<< "$cfgs")"
  fi
  if [[ ! -f extern/rnastructure_bridge/CMakeLists.txt ]]; then
    echo "Testing without RNAstructure because extern/rnastructure_bridge is missing"
    cfgs="$(sed -E 's/(^| )--rnastructure( |$)/\1--no-rnastructure\2/' <<< "$cfgs")"
  fi
  # With an option turned off, a configuration can be in the list more than once.
  # Compare sorted options to remove the repeats, since the options of the extra
  # builds above are in another order.
  cfgs="$(while read -r cfg; do
    printf '%s\t%s\n' "$(tr ' ' '\n' <<< "$cfg" | sort | paste -sd ' ')" "$cfg"
  done <<< "$cfgs" | awk -F '\t' '!seen[$1]++ { print $2 }')"
  # Suspend work when there is not enough memory for another job.
  printf '%s\n' "$cfgs" | parallel --halt soon,fail=1 --jobs "$test_jobs" \
    --memsuspend "$test_job_memory" just _test-build {}

_test-build *args:
  #!/usr/bin/env bash
  set -euo pipefail
  read -r -a cfg_args <<< "$*"
  # Build configurations run in parallel. Keep each build serial by
  # default so the two levels do not multiply the host CPU and memory use.
  export MRNA_BUILD_JOBS="${MRNA_BUILD_JOBS:-1}"
  prefix_args=()
  if [[ "${MRNA_CLEAN_BUILDS:-0}" == 1 ]]; then
    prefix="$(mktemp -d)"
    trap 'rm -rf "$prefix"' EXIT
    prefix_args=(--prefix "$prefix")
  fi
  /usr/bin/time --format="$* | %e s | %M KiB max RSS" \
    uv run python -m rnapy.run build --test \
      "${prefix_args[@]}" "${cfg_args[@]}" > /dev/null

docker-test:
  docker build --tag memerna:ubuntu-min .
  docker run --rm memerna:ubuntu-min

bench:
  # Run benchmarks.
  uv run python -m rnapy.run build --bench --bench-output \
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
    # With random models, the code for T12 is the same as for T04, and only the
    # parameter files differ. Fuzzing with --random-seeds does not load them.
    energy_models=(t04 t22)
    ctds=(none d2 no-coax all)
    random_model_ranges=(
      "--random-min-energy -10.0 --random-max-energy 10.0"
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
          # T22 is only supported by the stack backend, which does not support
          # D2. The fuzzer would exit successfully without testing anything.
          if [[ "$energy_model" == t22 && "$ctd" == d2 ]]; then continue; fi
          # The stack backend has no partition function for T22.
          if [[ "$energy_model" == t22 && "$base_cfg" == *--pfn* ]]; then continue; fi
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
  uv run python -m rnapy.run afl-fuzz-min --compiler=afl-fast --kind=relwithdebinfo \
    --mfe --mfe-table --subopt "$@" "${crashes[@]}"

afl-fuzz *args: afl-setup
  uv run python -m rnapy.run afl-fuzz --compiler=afl-fast --kind=relwithdebinfo \
    --mfe --mfe-table --subopt --random-pf "$@"

afl-run *args:
  uv run python -m rnapy.run afl-fuzz-run --compiler=afl-fast --kind=relwithdebinfo \
    --mfe --mfe-table --subopt --random-pf "$@"

fix:
  pre-commit run --all-files

update:
  uvx uv-upx upgrade run
  uv lock --upgrade

check:
  uv lock --check
  uv sync --locked
  uv run --locked uv pip check
  uv run --locked pytest
  uv run --locked ty check rnapy rnapy_tests
  uv run --locked ruff check rnapy rnapy_tests
  uv run --locked ruff format --check rnapy rnapy_tests
  pre-commit run --all-files
