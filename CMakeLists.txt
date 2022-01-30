cmake_minimum_required(VERSION 3.12)

project(memerna VERSION 0.1 LANGUAGES CXX)

# Options
option(USE_MPFR "Use MPFR for higher precision floats" OFF)
option(FLOAT_BITS "Number of bits to use for floats.
    Valid values are 32, 64, 80, or any value if MPFR is enabled." 64)
option(USE_RNASTRUCTURE "Whether to build in RNAstructure" OFF)
option(USE_IWYU "Whether to build with include-what-you-use" OFF)

# Always have RNAstructure compiled as release.
if (USE_RNASTRUCTURE)
  set(SAVED_BUILD_TYPE ${CMAKE_BUILD_TYPE})
  set(CMAKE_BUILD_TYPE Release)
  add_subdirectory("extern/miles_rnastructure")
  set(CMAKE_BUILD_TYPE ${SAVED_BUILD_TYPE})
endif ()

# Flags for each build configuration.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -march=native -O3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O3")

# Use include-what-you-use if it exists.
if (USE_IWYU)
  find_program(IWYU NAMES include-what-you-use iwyu REQUIRED)
  set(IWYU ${IWYU}
    -Xiwyu --max_line_length=100
    -Xiwyu --no_fwd_decls
    -Xiwyu --cxx17ns)
endif ()

# Include things to find libraries.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/sanitizers-cmake/cmake)
find_package(Sanitizers)
find_package(Threads)
find_package(MPFR)
find_package(Boost)

# Fail loudly.
# Note we don't support msan because it requires an instrumented libc++.
# See https://github.com/google/sanitizers/wiki/MemorySanitizerLibcxxHowTo
if (SANITIZE_ADDRESS AND NOT ASan_FLAG_DETECTED)
  message(FATAL_ERROR "requested asan build but not available on current compiler")
endif ()
if (SANITIZE_THREAD AND NOT TSan_FLAG_DETECTED)
  message(FATAL_ERROR "requested tsan build but not available on current compiler")
endif ()
if (SANITIZE_UNDEFINED AND NOT UBSan_FLAG_DETECTED)
  message(FATAL_ERROR "requested ubsan build but not available on current compiler")
endif ()

# Check for MPFR build.
if (USE_MPFR)
  if (NOT MPFR_FOUND)
    message(FATAL_ERROR "mpfr builds require mpfr")
  endif ()
  if (NOT Boost_FOUND)
    message(FATAL_ERROR "mpfr builds require boost")
  endif ()
  message(STATUS "Building using MPFR for floats.")
endif ()

# Enable tests.
set(THREADS_PREFER_PTHREAD_FLAG ON)  # Gtest needs this.
# Download and unpack googletest at configure time
configure_file(cmake/DownloadGTest.cmake googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
if (result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif ()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
if (result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif ()

add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
    ${CMAKE_CURRENT_BINARY_DIR}/googletest-build EXCLUDE_FROM_ALL)
include(GoogleTest)
enable_testing()

# Apply compile definitions here so they only apply to memerna code.
add_compile_options(-Wall -Wextra -Wconversion -Wno-sign-conversion -pedantic)

# Enable stacktraces if we have boost.
find_library(LIBBACKTRACE NAMES backtrace libbacktrace)
if (LIBBACKTRACE)
  message(STATUS "Using Boost stacktrace with libbacktrace - ${LIBBACKTRACE}")
  add_compile_definitions(BOOST_STACKTRACE_USE_BACKTRACE)
endif ()

# Source set definitions.
set(MEMERNA_SOURCE
    src/bridge/bridge.cpp src/bridge/bridge.h
    src/bridge/memerna.cpp src/bridge/memerna.h
    src/compute/boltz_dp.h
    src/compute/brute/brute.cpp src/compute/brute/brute.h
    src/compute/brute/alg.cpp src/compute/brute/alg.h
    src/compute/brute/config.cpp src/compute/brute/config.h
    src/compute/dp.h
    src/compute/energy/boltzmann_model.cpp src/compute/energy/boltzmann_model.h
    src/compute/energy/boltzmann_precomp.cpp src/compute/energy/boltzmann_precomp.h
    src/compute/energy/branch.cpp src/compute/energy/branch.h
    src/compute/energy/config.cpp src/compute/energy/config.h
    src/compute/energy/energy.cpp src/compute/energy/energy.h
    src/compute/energy/model.cpp src/compute/energy/model.h
    src/compute/energy/precomp.cpp src/compute/energy/precomp.h
    src/compute/energy/structure.cpp src/compute/energy/structure.h
    src/compute/mfe/mfe.cpp src/compute/mfe/mfe.h
    src/compute/mfe/mfe0.cpp
    src/compute/mfe/mfe1.cpp
    src/compute/mfe/mfe2.cpp
    src/compute/mfe/mfe3.cpp
    src/compute/partition/exterior.cpp
    src/compute/partition/partition.cpp src/compute/partition/partition.h
    src/compute/partition/partition0.cpp
    src/compute/partition/partition1.cpp
    src/compute/subopt/config.cpp src/compute/subopt/config.h
    src/compute/subopt/subopt.h
    src/compute/subopt/subopt0.cpp src/compute/subopt/subopt0.h
    src/compute/subopt/subopt1.cpp src/compute/subopt/subopt1.h
    src/compute/traceback/traceback.cpp src/compute/traceback/traceback.h
    src/ctx/config.cpp src/ctx/config.h
    src/ctx/ctx.cpp src/ctx/ctx.h
    src/fuzz/config.cpp src/fuzz/config.h
    src/fuzz/fuzzer.cpp src/fuzz/fuzzer.h
    src/model/base.cpp src/model/base.h
    src/model/ctd.cpp src/model/ctd.h
    src/model/model.h
    src/model/primary.cpp src/model/primary.h
    src/model/secondary.cpp src/model/secondary.h
    src/util/argparse.cpp src/util/argparse.h
    src/programs/print.cpp src/programs/print.h
    src/util/array.h
    src/util/splaymap.h
    src/util/string.cpp src/util/string.h
)
set(TEST_SOURCE
    tests/common_test.h
    tests/compute/brute/brute_test.cpp
    tests/compute/energy/energy_test.cpp
    tests/compute/mfe/mfe_test.cpp
    tests/model/base_test.cpp
    tests/model/ctds_test.cpp
    tests/model/parsing_test.cpp
    tests/util/splaymap_test.cpp
)

# Use RNAstructure if necessary
if (USE_RNASTRUCTURE)
  set(MEMERNA_SOURCE ${MEMERNA_SOURCE}
      src/bridge/rnastructure.cpp src/bridge/rnastructure.h
  )
endif ()

add_library(memerna ${MEMERNA_SOURCE})

add_executable(brute src/programs/brute.cpp)
add_executable(efn src/programs/efn.cpp)
add_executable(fold src/programs/fold.cpp)
add_executable(fuzz src/programs/fuzz.cpp)
add_executable(harness src/programs/harness.cpp)
add_executable(partition src/programs/partition.cpp)
add_executable(run_tests ${TEST_SOURCE} tests/programs/run_tests.cpp)
add_executable(splay_explorer tests/programs/splay_explorer.cpp)
add_executable(subopt src/programs/subopt.cpp)

target_include_directories(memerna PUBLIC src)
target_include_directories(run_tests PRIVATE tests)

target_link_libraries(brute PRIVATE memerna)
target_link_libraries(efn PRIVATE memerna)
target_link_libraries(fold PRIVATE memerna)
target_link_libraries(fuzz PRIVATE memerna)
target_link_libraries(harness PRIVATE memerna)
target_link_libraries(partition PRIVATE memerna)
target_link_libraries(run_tests PRIVATE memerna Threads::Threads gtest)
target_link_libraries(splay_explorer PRIVATE memerna)
target_link_libraries(subopt PRIVATE memerna)

# Use Boost if available
if (Boost_FOUND)
  target_include_directories(memerna PUBLIC ${Boost_INCLUDE_DIRS})
  target_compile_definitions(memerna PUBLIC USE_BOOST)
  target_link_libraries(memerna PUBLIC ${CMAKE_DL_LIBS})  # For Boost stacktrace

  if (LIBBACKTRACE)
    target_link_libraries(memerna PUBLIC ${LIBBACKTRACE})
  endif ()
endif ()

# Use MPFR if necessary
if (USE_MPFR)
  target_include_directories(memerna PUBLIC ${MPFR_INCLUDE_DIRS})
  target_compile_definitions(memerna PUBLIC USE_MPFR)
  target_link_libraries(memerna PUBLIC ${MPFR_LIBRARIES})
endif ()

# Use RNAstructure if necessary
if (USE_RNASTRUCTURE)
  target_include_directories(memerna SYSTEM PUBLIC extern) # Suppress warnings.
  target_compile_definitions(memerna PUBLIC USE_RNASTRUCTURE)
  target_link_libraries(memerna PUBLIC miles_rnastructure)
endif ()

# Add sanitizers
add_sanitizers(memerna)
add_sanitizers(brute)
add_sanitizers(efn)
add_sanitizers(fold)
add_sanitizers(fuzz)
add_sanitizers(harness)
add_sanitizers(partition)
add_sanitizers(run_tests)
add_sanitizers(splay_explorer)
add_sanitizers(subopt)

# Add iwyu to targets if requested.
if (USE_IWYU)
  set_target_properties(memerna PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(brute PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(efn PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(fold PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(fuzz PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(harness PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(partition PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(run_tests PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(splay_explorer PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
  set_target_properties(subopt PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "${IWYU}")
endif ()
