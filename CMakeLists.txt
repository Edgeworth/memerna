cmake_minimum_required(VERSION 3.12)

project(memerna VERSION 0.1 LANGUAGES CXX)

# Supress warning messages about unused variables:
set(SUPPRESS_UNUSED "${CMAKE_C_COMPILER}")

# Options
option(USE_MPFR "Use MPFR for higher precision floats" OFF)
option(USE_RNASTRUCTURE "Whether to build in RNAstructure" ON)
option(USE_IWYU "Whether to build with include-what-you-use" OFF)
option(USE_LTO "Whether to build with LTO" OFF)
set(FLOAT_PRECISION 15 CACHE STRING "Number of significant digits to use for floats.
    Valid values are 6 (float32), 15 (float64), 18 (float80), or any value if MPFR is enabled.")
set(ENERGY_PRECISION 2 CACHE STRING
  "Number of decimal places to use for energy calculations")

# Always have RNAstructure compiled as release.
if(USE_RNASTRUCTURE)
  set(SAVED_BUILD_TYPE ${CMAKE_BUILD_TYPE})
  set(CMAKE_BUILD_TYPE Release)
  add_subdirectory("extern/rnastructure_bridge")
  set(CMAKE_BUILD_TYPE ${SAVED_BUILD_TYPE})
endif()

# Flags for each build configuration.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -march=native -O3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O3")

# Enable LTO
if(USE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported()
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# Use include-what-you-use if it exists.
if(USE_IWYU)
  find_program(IWYU NAMES include-what-you-use iwyu REQUIRED)
  set(IWYU ${IWYU}
    -Xiwyu --max_line_length=100
    -Xiwyu --no_fwd_decls
    -Xiwyu --cxx17ns)
endif()

# Include things to find libraries.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/sanitizers-cmake/cmake)
include(get_cpm) # get cmake package manager

# Find packages
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

find_package(Threads REQUIRED)
find_package(Boost REQUIRED)
find_package(Sanitizers)
find_package(MPFR)

# Enable stacktraces if we have boost.
find_library(LIBBACKTRACE NAMES backtrace libbacktrace)

if(LIBBACKTRACE)
  message(STATUS "Using Boost stacktrace with libbacktrace - ${LIBBACKTRACE}")
  add_compile_definitions(BOOST_STACKTRACE_USE_BACKTRACE)
endif()

# Get rarer packages using CPM.

# Google benchmark
CPMAddPackage(
  NAME benchmark
  GITHUB_REPOSITORY google/benchmark
  VERSION 1.9.1
  OPTIONS "BENCHMARK_ENABLE_TESTING Off"
)

# fmt
CPMAddPackage("gh:fmtlib/fmt#11.0.2")

# spdlog should reuse existing fmtlib.
CPMAddPackage(
  GITHUB_REPOSITORY gabime/spdlog
  VERSION 1.15.0
  OPTIONS "SPDLOG_FMT_EXTERNAL 1"
)

# gtest
CPMAddPackage(
  NAME googletest
  GITHUB_REPOSITORY google/googletest
  VERSION 1.15.2
  OPTIONS "INSTALL_GTEST OFF" "gtest_force_shared_crt"
)

# Fail loudly.
# Note we don't support msan because it requires an instrumented libc++.
# See https://github.com/google/sanitizers/wiki/MemorySanitizerLibcxxHowTo
if(SANITIZE_ADDRESS AND NOT ASan_FLAG_DETECTED)
  message(FATAL_ERROR "requested asan build but not available on current compiler")
endif()

if(SANITIZE_THREAD AND NOT TSan_FLAG_DETECTED)
  message(FATAL_ERROR "requested tsan build but not available on current compiler")
endif()

if(SANITIZE_UNDEFINED AND NOT UBSan_FLAG_DETECTED)
  message(FATAL_ERROR "requested ubsan build but not available on current compiler")
endif()

# Check for MPFR build.
if(USE_MPFR)
  if(NOT MPFR_FOUND)
    message(FATAL_ERROR "mpfr builds require mpfr")
  endif()

  if(NOT Boost_FOUND)
    message(FATAL_ERROR "mpfr builds require boost")
  endif()

  message(STATUS "Building using MPFR for floats.")
endif()

# Enable tests.
set(THREADS_PREFER_PTHREAD_FLAG ON) # Gtest needs this.
enable_testing()

# Apply compile definitions here so they only apply to memerna code.
add_compile_options(-Wall -Wextra -Wconversion -Wno-sign-conversion -pedantic)

# Source set definitions.
set(MEMERNA_SOURCE
  src/api/bridge/bridge.cpp
  src/api/bridge/bridge.h
  src/api/bridge/memerna.cpp
  src/api/bridge/memerna.h
  src/api/brute/brute_cfg.cpp
  src/api/brute/brute_cfg.h
  src/api/ctx/backend_cfg.cpp
  src/api/ctx/backend_cfg.h
  src/api/ctx/backend.cpp
  src/api/ctx/backend.h
  src/api/ctx/ctx_cfg.cpp
  src/api/ctx/ctx_cfg.h
  src/api/ctx/ctx.cpp
  src/api/ctx/ctx.h
  src/api/energy/energy_cfg.cpp
  src/api/energy/energy_cfg.h
  src/api/energy/energy.cpp
  src/api/energy/energy.h
  src/api/mfe.h
  src/api/options.h
  src/api/pfn.h
  src/api/subopt/subopt_cfg.cpp
  src/api/subopt/subopt_cfg.h
  src/api/subopt/subopt.h
  src/api/trace/trace_cfg.cpp
  src/api/trace/trace_cfg.h
  src/api/trace/trace.h
  src/backends/base/energy/boltz_model.h
  src/backends/base/energy/boltz_precomp.cpp
  src/backends/base/energy/boltz_precomp.h
  src/backends/base/energy/model.cpp
  src/backends/base/energy/model.h
  src/backends/base/energy/precomp.cpp
  src/backends/base/energy/precomp.h
  src/backends/base/mfe/mfe_debug.cpp
  src/backends/base/mfe/mfe_exterior.cpp
  src/backends/base/mfe/mfe_lyngso_sparse_opt.cpp
  src/backends/base/mfe/mfe_opt.cpp
  src/backends/base/mfe/mfe_sparse_opt.cpp
  src/backends/base/mfe/mfe.h
  src/backends/base/pfn/pfn_debug.cpp
  src/backends/base/pfn/pfn_exterior.cpp
  src/backends/base/pfn/pfn_opt.cpp
  src/backends/base/pfn/pfn.h
  src/backends/base/subopt/subopt_debug.cpp
  src/backends/base/subopt/subopt_debug.h
  src/backends/base/subopt/subopt_iterative.cpp
  src/backends/base/subopt/subopt_iterative.h
  src/backends/base/subopt/subopt_persistent.cpp
  src/backends/base/subopt/subopt_persistent.h
  src/backends/base/trace/trace.cpp
  src/backends/base/trace/trace.h
  src/backends/baseopt/energy/boltz_model.h
  src/backends/baseopt/energy/boltz_precomp.cpp
  src/backends/baseopt/energy/boltz_precomp.h
  src/backends/baseopt/energy/model.cpp
  src/backends/baseopt/energy/model.h
  src/backends/baseopt/energy/precomp.cpp
  src/backends/baseopt/energy/precomp.h
  src/backends/baseopt/mfe/mfe_debug.cpp
  src/backends/baseopt/mfe/mfe_exterior.cpp
  src/backends/baseopt/mfe/mfe_lyngso_sparse_opt.cpp
  src/backends/baseopt/mfe/mfe_opt.cpp
  src/backends/baseopt/mfe/mfe_sparse_opt.cpp
  src/backends/baseopt/mfe/mfe.h
  src/backends/baseopt/pfn/pfn_debug.cpp
  src/backends/baseopt/pfn/pfn_exterior.cpp
  src/backends/baseopt/pfn/pfn_opt.cpp
  src/backends/baseopt/pfn/pfn.h
  src/backends/baseopt/subopt/subopt_debug.cpp
  src/backends/baseopt/subopt/subopt_debug.h
  src/backends/baseopt/subopt/subopt_iterative.cpp
  src/backends/baseopt/subopt/subopt_iterative.h
  src/backends/baseopt/subopt/subopt_persistent.cpp
  src/backends/baseopt/subopt/subopt_persistent.h
  src/backends/baseopt/trace/trace.cpp
  src/backends/baseopt/trace/trace.h
  src/backends/brute/alg.cpp
  src/backends/brute/alg.h
  src/backends/brute/brute.cpp
  src/backends/brute/brute.h
  src/backends/common/base/boltz_model_base.h
  src/backends/common/base/boltz_precomp_base.h
  src/backends/common/base/branch.cpp
  src/backends/common/base/branch.h
  src/backends/common/base/dp.h
  src/backends/common/base/model_base.h
  src/backends/common/base/parse.h
  src/backends/common/base/precomp_base.h
  src/backends/common/base/pseudofree_model.cpp
  src/backends/common/base/pseudofree_model.h
  src/backends/common/boltz.h
  src/backends/common/branch.cpp
  src/backends/common/branch.h
  src/backends/common/energy.h
  src/backends/common/model_mixin.h
  src/backends/common/parse.cpp
  src/backends/common/parse.h
  src/backends/stack/energy/boltz_model.cpp
  src/backends/stack/energy/boltz_model.h
  src/backends/stack/energy/model.cpp
  src/backends/stack/energy/model.h
  src/backends/stack/mfe/mfe_debug.cpp
  src/backends/stack/mfe/mfe_exterior.cpp
  src/backends/stack/mfe/mfe.h
  src/backends/stack/subopt/subopt_iterative.cpp
  src/backends/stack/subopt/subopt_iterative.h
  src/backends/stack/subopt/subopt_persistent.cpp
  src/backends/stack/subopt/subopt_persistent.h
  src/backends/stack/trace/trace.cpp
  src/backends/stack/trace/trace.h
  src/programs/print.cpp
  src/programs/print.h
  src/fuzz/fuzz_cfg.cpp
  src/fuzz/fuzz_cfg.h
  src/fuzz/fuzz_harness.cpp
  src/fuzz/fuzz_harness.h
  src/fuzz/fuzz_invocation.cpp
  src/fuzz/fuzz_invocation.h
  src/model/base.cpp
  src/model/base.h
  src/model/constants.h
  src/model/ctd.cpp
  src/model/ctd.h
  src/model/energy.cpp
  src/model/energy.h
  src/model/pfn.h
  src/model/primary.cpp
  src/model/primary.h
  src/model/secondary.cpp
  src/model/secondary.h
  src/model/structure.cpp
  src/model/structure.h
  src/util/argparse.cpp
  src/util/argparse.h
  src/util/array.h
  src/util/container.h
  src/util/error.cpp
  src/util/error.h
  src/util/float.h
  src/util/splaymap.h
  src/util/string.cpp
  src/util/string.h
  src/util/util.h
)
set(TEST_SOURCE
  src/tests/backends/energy_test.cpp
  src/tests/backends/baseopt/energy_test.cpp
  src/tests/backends/stack/energy_test.cpp
  src/tests/backends/base/energy_test.cpp
  src/tests/model/base_test.cpp
  src/tests/model/parsing_test.cpp
  src/tests/model/model_test.cpp
  src/tests/models/t22/energy_test.cpp
  src/tests/models/t22/mfe_test.cpp
  src/tests/models/t12/energy_test.cpp
  src/tests/models/t12/subopt_test.cpp
  src/tests/models/t12/pfn_test.cpp
  src/tests/models/t12/mfe_test.cpp
  src/tests/models/t04/energy_test.cpp
  src/tests/models/t04/subopt_test.cpp
  src/tests/models/t04/pfn_test.cpp
  src/tests/models/t04/mfe_test.cpp
  src/tests/util/splaymap_test.cpp
  src/tests/util.cpp
  src/tests/util.h
  src/tests/init.cpp
  src/tests/init.h
)
set(BENCHMARK_SOURCE
  src/tests/init.cpp
  src/tests/init.h
)

# Use RNAstructure if necessary
if(USE_RNASTRUCTURE)
  set(MEMERNA_SOURCE ${MEMERNA_SOURCE}
    src/api/bridge/rnastructure.cpp src/api/bridge/rnastructure.h
  )
endif()

add_library(memerna ${MEMERNA_SOURCE})

add_executable(brute src/programs/brute.cpp)
add_executable(efn src/programs/efn.cpp)
add_executable(fold src/programs/fold.cpp)
add_executable(fuzz src/programs/fuzz/fuzz.cpp)
add_executable(fuzz_afl src/programs/fuzz/fuzz_afl.cpp)
add_executable(harness src/programs/harness.cpp)
add_executable(partition src/programs/partition.cpp)
add_executable(run_benchmark ${BENCHMARK_SOURCE} src/tests/programs/run_benchmark.cpp)
add_executable(run_tests ${TEST_SOURCE} src/tests/programs/run_tests.cpp)
add_executable(splay_explorer src/tests/programs/splay_explorer.cpp)
add_executable(subopt src/programs/subopt.cpp)

target_include_directories(memerna PUBLIC src)

target_link_libraries(memerna PUBLIC fmt::fmt spdlog::spdlog)
target_link_libraries(brute PRIVATE memerna)
target_link_libraries(efn PRIVATE memerna)
target_link_libraries(fold PRIVATE memerna)
target_link_libraries(fuzz PRIVATE memerna)
target_link_libraries(fuzz_afl PRIVATE memerna)
target_link_libraries(harness PRIVATE memerna)
target_link_libraries(partition PRIVATE memerna)
target_link_libraries(run_benchmark PRIVATE memerna Threads::Threads benchmark::benchmark)
target_link_libraries(run_tests PRIVATE memerna Threads::Threads GTest::gtest)
target_link_libraries(splay_explorer PRIVATE memerna)
target_link_libraries(subopt PRIVATE memerna)

# Use Boost if available
if(Boost_FOUND)
  target_include_directories(memerna PUBLIC ${Boost_INCLUDE_DIRS})
  target_compile_definitions(memerna PUBLIC USE_BOOST)
  target_link_libraries(memerna PUBLIC ${CMAKE_DL_LIBS}) # For Boost stacktrace

  if(LIBBACKTRACE)
    target_link_libraries(memerna PUBLIC ${LIBBACKTRACE})
  endif()
endif()

# Use MPFR if necessary
if(USE_MPFR)
  target_include_directories(memerna PUBLIC ${MPFR_INCLUDE_DIRS})
  target_compile_definitions(memerna PUBLIC USE_MPFR)
  target_link_libraries(memerna PUBLIC ${MPFR_LIBRARIES})
endif()

# Use RNAstructure if necessary
if(USE_RNASTRUCTURE)
  target_include_directories(memerna SYSTEM PUBLIC extern) # Suppress warnings.
  target_compile_definitions(memerna PUBLIC USE_RNASTRUCTURE)
  target_link_libraries(memerna PUBLIC rnastructure_bridge)
endif()

# Set float precision
target_compile_definitions(memerna PUBLIC FLOAT_PRECISION=${FLOAT_PRECISION})

# Set energy precision
target_compile_definitions(memerna PUBLIC ENERGY_PRECISION=${ENERGY_PRECISION})

# Add sanitizers
add_sanitizers(memerna)
add_sanitizers(brute)
add_sanitizers(efn)
add_sanitizers(fold)
add_sanitizers(fuzz)
add_sanitizers(fuzz_afl)
add_sanitizers(harness)
add_sanitizers(partition)
add_sanitizers(run_benchmark)
add_sanitizers(run_tests)
add_sanitizers(splay_explorer)
add_sanitizers(subopt)

# Add iwyu to targets if requested.
if(USE_IWYU)
  set_target_properties(memerna PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(brute PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(efn PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(fold PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(fuzz PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(fuzz_afl PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(harness PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(partition PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(run_benchmark PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(run_tests PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(splay_explorer PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(subopt PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
endif()

# Copy data assets to build directory, if the build succeeds:
add_custom_command(TARGET memerna POST_BUILD COMMAND ${CMAKE_COMMAND}
  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data/
  $<TARGET_FILE_DIR:memerna>/data)
