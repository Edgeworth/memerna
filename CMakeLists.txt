cmake_minimum_required(VERSION 3.2)

project(memerna)

set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O2")
file(GLOB MEMERNA_SOURCE "src/*.cpp" "src/*.h")
file(GLOB TEST_SOURCE "tests/*.cpp" "tests/*.h")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FindUBSan)
include(FindASan)
include(FindMSan)
include(FindThreads)

include_directories(src)

add_library(memerna ${MEMERNA_SOURCE})

add_executable(compute_energy src/programs/compute_energy.cpp)
add_executable(batch_energy src/programs/batch_energy.cpp)
add_executable(fold src/programs/fold.cpp)
add_executable(brute_fold src/programs/brute_fold.cpp)
add_executable(fuzz src/programs/fuzz.cpp)

add_subdirectory(lib/gtest)
enable_testing()
add_executable(run_tests ${TEST_SOURCE} ${MEMERNA_SOURCE})

target_link_libraries(run_tests memerna Threads::Threads gtest)
target_link_libraries(compute_energy memerna)
target_link_libraries(batch_energy memerna)
target_link_libraries(fold memerna)
target_link_libraries(brute_fold memerna)
target_link_libraries(fuzz memerna)
