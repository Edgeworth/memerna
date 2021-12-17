cmake_minimum_required(VERSION 3.12)

project(memerna VERSION 0.1 LANGUAGES CXX)

# Options
option(PARTITION_MPFR "Use higher precision with MPFR for the partition function" OFF)

# Always have RNAstructure compiled as release.
set(SAVED_BUILD_TYPE ${CMAKE_BUILD_TYPE})
set(CMAKE_BUILD_TYPE Release)
add_subdirectory("extern/miles_rnastructure")
set(CMAKE_BUILD_TYPE ${SAVED_BUILD_TYPE})

# Flags for each build configuration.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_ASAN "${CMAKE_CXX_FLAGS_ASAN} -march=native -O3")
set(CMAKE_CXX_FLAGS_UBSAN "${CMAKE_CXX_FLAGS_UBSAN} -march=native -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -march=native -O3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O3")

# Include things to find libraries.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
find_package(UBSan)
find_package(ASan)
find_package(Threads)
find_package(MPFR)
find_package(Boost)

# Fail loudly.
string(TOLOWER "${CMAKE_BUILD_TYPE}" SAVED_BUILD_TYPE)
if ("${SAVED_BUILD_TYPE}" STREQUAL "asan" AND NOT HAVE_ADDRESS_SANITIZER)
  message(FATAL_ERROR "requested asan build but not available on current compiler")
endif ()
if ("${SAVED_BUILD_TYPE}" STREQUAL "ubsan" AND NOT HAVE_UNDEFINED_BEHAVIOR_SANITIZER)
  message(FATAL_ERROR "requested ubsan build but not available on current compiler")
endif ()

# Check for MPFR build.
if (PARTITION_MPFR)
  if (NOT MPFR_FOUND)
    message(FATAL_ERROR "mpfr builds require mpfr")
  endif ()
  if (NOT Boost_FOUND)
    message(FATAL_ERROR "mpfr builds require boost")
  endif ()
  include_directories(${Boost_INCLUDE_DIRS} ${MPFR_INCLUDE_DIRS})
  add_definitions(-DPARTITION_MPFR)
  message("Building using MPFR for the partition function.")
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

# Source set definitions.
set(MEMERNA_SOURCE
    src/common.cpp src/common.h
    src/compute/constants.h
    src/compute/dp.h
    src/compute/energy/energy.cpp src/compute/energy/energy.h
    src/compute/energy/fast_energy.cpp src/compute/energy/fast_energy.h
    src/compute/energy/globals.cpp src/compute/energy/globals.h
    src/compute/energy/internal.cpp src/compute/energy/internal.h
    src/compute/energy/load_model.cpp src/compute/energy/load_model.h
    src/compute/energy/model.cpp src/compute/energy/model.h
    src/compute/energy/structure.cpp src/compute/energy/structure.h
    src/compute/mfe/globals.cpp src/compute/mfe/globals.h
    src/compute/mfe/brute.cpp src/compute/mfe/brute.h
    src/compute/mfe/mfe.cpp src/compute/mfe/mfe.h
    src/compute/mfe/mfe0.cpp
    src/compute/mfe/mfe1.cpp
    src/compute/mfe/mfe2.cpp
    src/compute/mfe/mfe3.cpp
    src/compute/partition/exterior.cpp
    src/compute/partition/globals.cpp src/compute/partition/globals.h
    src/compute/partition/internal.cpp
    src/compute/partition/brute.cpp src/compute/partition/brute.h
    src/compute/partition/partition.h
    src/compute/partition/partition0.cpp
    src/compute/partition/partition1.cpp
    src/compute/partition/probability.cpp
    src/compute/subopt/brute.cpp src/compute/subopt/brute.h
    src/compute/subopt/subopt.h
    src/compute/subopt/subopt0.cpp src/compute/subopt/subopt0.h
    src/compute/subopt/subopt1.cpp src/compute/subopt/subopt1.h
    src/compute/traceback/traceback.cpp src/compute/traceback/traceback.h
    src/compute/brute.cpp src/compute/brute.h
    src/model/base.cpp src/model/base.h
    src/model/context.cpp src/model/context.h
    src/model/gen.cpp src/model/gen.h
    src/model/globals.cpp src/model/globals.h
    src/model/parsing.cpp src/model/parsing.h
    src/model/structure.h
    src/util/argparse.cpp src/util/argparse.h
    src/util/array.h
    src/util/splaymap.h
    src/util/string.cpp src/util/string.h
)
set(BRIDGE_SOURCE
    src/bridge/bridge.cpp src/bridge/bridge.h
    src/bridge/memerna.cpp src/bridge/memerna.h
    src/bridge/rnastructure.cpp src/bridge/rnastructure.h
)
set(TEST_SOURCE
    tests/base_test.cpp
    tests/brute_test.cpp
    tests/common_test.cpp
    tests/common_test.h
    tests/ctds_test.cpp
    tests/energy_test.cpp
    tests/mfe_test.cpp
    tests/parsing_test.cpp
    tests/splaymap_test.cpp
)

include_directories(src)
# Suppress warnings.
include_directories(SYSTEM extern)

# Libraries.
add_library(memerna ${MEMERNA_SOURCE})
add_library(bridge ${BRIDGE_SOURCE})

# Link MPFR if necessary
if (PARTITION_MPFR)
  target_link_libraries(memerna ${MPFR_LIBRARIES})
endif ()

add_executable(efn src/programs/efn.cpp)
add_executable(fold src/programs/fold.cpp)
add_executable(subopt src/programs/subopt.cpp)
add_executable(fuzz src/programs/fuzz.cpp)
add_executable(harness src/programs/harness.cpp)
add_executable(partition src/programs/partition.cpp)
add_executable(run_tests ${TEST_SOURCE} tests/programs/run_tests.cpp)
add_executable(splay_explorer tests/programs/splay_explorer.cpp)

target_include_directories(run_tests PRIVATE tests)

target_link_libraries(efn memerna)
target_link_libraries(fold memerna)
target_link_libraries(subopt memerna)
target_link_libraries(fuzz bridge memerna miles_rnastructure)
target_link_libraries(harness bridge memerna miles_rnastructure)
target_link_libraries(partition bridge memerna miles_rnastructure)
target_link_libraries(run_tests memerna Threads::Threads gtest)
target_link_libraries(splay_explorer memerna)
