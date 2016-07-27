cmake_minimum_required(VERSION 3.2)

project(memerna)

set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O2")
file(GLOB MEMERNA_SOURCE "src/*.cpp" "src/*.h")
list(REMOVE_ITEM MEMERNA_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
file(GLOB TEST_SOURCE "tests/*.cpp" "tests/*.h")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FindUBSan)
include(FindASan)
include(FindMSan)
include(FindThreads)

include_directories(src)

add_executable(memerna ${MEMERNA_SOURCE} src/main.cpp)

add_subdirectory(lib/gtest)
enable_testing()
add_executable(run_tests ${TEST_SOURCE} ${MEMERNA_SOURCE})
target_link_libraries(run_tests Threads::Threads gtest)
