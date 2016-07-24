cmake_minimum_required(VERSION 3.2)

project(memerna)

set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O3")
message(${CMAKE_CXX_FLAGS_RELEASE})
file(GLOB MEMERNA_SOURCE "src/*.cpp")
file(GLOB MEMERNA_HEADERS "src/*.h")
list(REMOVE_ITEM MEMERNA_SOURCE src/main.cpp)
file(GLOB TEST_SOURCE "tests/*.cpp")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FindUBSan)
include(FindASan)
include(FindMSan)
include(FindThreads)

include_directories(src)

add_library(memerna_lib ${MEMERNA_SOURCE} ${MEMERNA_HEADERS} src/common.cpp)
add_executable(memerna src/main.cpp)
target_link_libraries(memerna memerna_lib)

add_subdirectory(lib/gtest)
enable_testing()
add_executable(run_tests ${TEST_SOURCE} ${MEMERNA_HEADERS})
target_link_libraries(run_tests Threads::Threads gtest memerna_lib)
