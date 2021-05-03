cmake_minimum_required(VERSION 3.2)

project(memerna)

# Options
option(PARTITION_MPFR "Use higher precision with MPFR for the partition function" OFF)

# I want LTO to work if I decide to use it.
set(BUILD_SHARED_LIBS OFF)

# Always have RNAstructure compiled as release.
set(SAVED_BUILD_TYPE ${CMAKE_BUILD_TYPE})
set(CMAKE_BUILD_TYPE Release)
add_subdirectory("extern/miles_rnastructure")
set(CMAKE_BUILD_TYPE ${SAVED_BUILD_TYPE})

# Include things to find libraries.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
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

# Gold linker.
check_c_compiler_flag("-fuse-ld=gold" HAVE_GOLD)
if (HAVE_GOLD)
  message("Using gold linker")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=gold")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fuse-ld=gold")
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

# Flags for each build configuration.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -Wshadow \
  -fno-rtti -Wconversion -Wfloat-equal -Wformat-security -Wno-sign-conversion -pedantic")
set(CMAKE_CXX_FLAGS_ASAN "${CMAKE_CXX_FLAGS_ASAN} -march=native -O3")
set(CMAKE_CXX_FLAGS_UBSAN "${CMAKE_CXX_FLAGS_UBSAN} -march=native -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -fno-omit-frame-pointer -march=native -O3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O3")

# Source set definitions.
file(GLOB MEMERNA_SOURCE "src/*.cpp" "src/*.h" "src/energy/*.cpp"
    "src/energy/*.h" "src/fold/*.cpp" "src/fold/*.h" "src/partition/*.h" "src/partition/*.cpp")
file(GLOB BRIDGE_SOURCE "src/bridge/*.cpp" "src/bridge/*.h")
file(GLOB TEST_SOURCE "tests/*.cpp" "tests/*.h")

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
