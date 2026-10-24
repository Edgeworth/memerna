alias t := test
alias f := fix
alias u := update

default:
  @just --list

test:
  #!/usr/bin/env bash
  set -euo pipefail
  test_job_memory="${MRNA_TEST_JOB_MEMORY:-4G}"
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
  # Check various build configurations and run tests. Good to run before push.
  # Uses :::+ to match mpfr-rnastructure and no-mpfr-no-rnastructure etc to
  # save some compilations.
  # Necessary to compile with multiple compilers, as some issues have only shown
  # themselves on a specific compiler.
  cfgs="$(parallel --keep-order echo ::: \
    --kind=debug --kind=relwithdebinfo ::: --sanitizer=asan \
    --sanitizer=ubsan ::: --float-precision=15 --float-precision=18 :::+ \
    --energy-precision=1 --energy-precision=2 :::+ \
    --rnastructure --no-rnastructure :::+ --mpfr --no-mpfr ::: \
    "${test_compilers[@]}")"
  # Also check one build without a sanitizer, since sanitizers change codegen.
  for test_compiler in "${test_compilers[@]}"; do
    cfgs+=$'\n'"--kind=debug --sanitizer=none --float-precision=15"
    cfgs+=" --energy-precision=1 --rnastructure --mpfr $test_compiler"
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
  cfgs="$(awk '!seen[$0]++' <<< "$cfgs")"
  # Suspend work when there is not enough memory for another job.
  printf '%s\n' "$cfgs" | parallel --halt soon,fail=1 --jobs "$test_jobs" \
    --memsuspend "$test_job_memory" just _test-build {}

_test-build *args:
  #!/usr/bin/env bash
  set -euo pipefail
  read -r -a cfg_args <<< "{{args}}"
  # Build configurations run in parallel. Keep each build serial by
  # default so the two levels do not multiply the host CPU and memory use.
  export MRNA_BUILD_JOBS="${MRNA_BUILD_JOBS:-1}"
  prefix_args=()
  if [[ "${MRNA_CLEAN_BUILDS:-0}" == 1 ]]; then
    prefix="$(mktemp -d)"
    trap 'rm -rf "$prefix"' EXIT
    prefix_args=(--prefix "$prefix")
  fi
  /usr/bin/time --format="{{args}} | %e s | %M KiB max RSS" \
    poetry run python -m rnapy.run build --test \
      "${prefix_args[@]}" "${cfg_args[@]}" > /dev/null

docker-test:
  docker build --tag memerna:ubuntu-min .
  docker run --rm memerna:ubuntu-min

bench:
  # Run benchmarks.
  poetry run python -m rnapy.run build --bench --bench-output \
    ./data/benchmark.json --kind=release

fuzz $fuzz_exec:
  #!/usr/bin/env -S parallel --shebang --ungroup --verbose
  $fuzz_exec --mfe --mfe-table --subopt --pfn --energy-model t04 1 30
  $fuzz_exec --mfe --mfe-table --subopt --pfn --energy-model t04 1 200
  $fuzz_exec --mfe --mfe-table --subopt --random-models 1 200
  $fuzz_exec --mfe --mfe-table --subopt --pfn --backends base,baseopt --random-models 1 30
  $fuzz_exec --mfe --mfe-table --subopt --pfn --backends base --random-models --ctd none 1 30
  $fuzz_exec --mfe --mfe-table --subopt --pfn --backends base --random-models --ctd d2 1 30
  $fuzz_exec --mfe --mfe-table --subopt --pfn --backends base --random-models --ctd no-coax 1 30
  $fuzz_exec --mfe --mfe-table --subopt --pfn --backends base --random-models --ctd all 1 30
  $fuzz_exec --mfe --mfe-table --subopt --pfn --backends base --random-models --random-pf --ctd none 1 30
  $fuzz_exec --mfe --mfe-table --subopt --pfn --backends base --random-models --random-pf --ctd d2 1 30
  $fuzz_exec --mfe --mfe-table --subopt --pfn --backends base --random-models --random-pf --ctd no-coax 1 30
  $fuzz_exec --mfe --mfe-table --subopt --pfn --backends base --random-models --random-pf --ctd all 1 30
  $fuzz_exec --mfe --mfe-table --subopt --backends base,stack --random-models --ctd none 1 30
  $fuzz_exec --mfe --mfe-table --subopt --backends base --random-models --ctd d2 1 30
  $fuzz_exec --mfe --mfe-table --subopt --backends base,stack --random-models --ctd no-coax 1 30
  $fuzz_exec --mfe --mfe-table --subopt --backends base,stack --random-models --ctd all 1 30
  $fuzz_exec --mfe --mfe-table --subopt --backends base,stack --random-models --ctd none 1 200
  $fuzz_exec --mfe --mfe-table --subopt --backends base --random-models --ctd d2 1 200
  $fuzz_exec --mfe --mfe-table --subopt --backends base,stack --random-models --ctd no-coax 1 200
  $fuzz_exec --mfe --mfe-table --subopt --backends base,stack --random-models --ctd all 1 200
  $fuzz_exec --mfe --mfe-table --subopt --backends base,stack --random-models --random-pf --ctd none 1 30
  $fuzz_exec --mfe --mfe-table --subopt --backends base --random-models --random-pf --ctd d2 1 30
  $fuzz_exec --mfe --mfe-table --subopt --backends base,stack --random-models --random-pf --ctd no-coax 1 30
  $fuzz_exec --mfe --mfe-table --subopt --backends base,stack --random-models --random-pf --ctd all 1 30
  $fuzz_exec --mfe --mfe-table --subopt --backends base,stack --random-models --random-pf --ctd none 1 200
  $fuzz_exec --mfe --mfe-table --subopt --backends base --random-models --random-pf --ctd d2 1 200
  $fuzz_exec --mfe --mfe-table --subopt --backends base,stack --random-models --random-pf --ctd no-coax 1 200
  $fuzz_exec --mfe --mfe-table --subopt --backends base,stack --random-models --random-pf --ctd all 1 200

fix:
  pre-commit run --all-files

update:
  poetry run poetry up --latest
  poetry update
  pre-commit autoupdate
