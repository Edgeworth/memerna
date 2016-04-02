cmake_minimum_required(VERSION 3.2)

project(memerna)

set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")
file(GLOB MEMERNA_SOURCE "src/*.cpp")
file(GLOB MEMERNA_HEADERS "src/*.h")
list(REMOVE_ITEM MEMERNA_SOURCE src/main.cpp)
file(GLOB TEST_SOURCE "tests/*.cpp")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FindUBSan)
include(FindASan)
include(FindMSan)
include(FindGTest)
include(FindThreads)

add_library(memerna_lib ${MEMERNA_SOURCE} ${MEMERNA_HEADERS})
add_executable(memerna src/main.cpp)
target_link_libraries(memerna memerna_lib)

enable_testing()
add_executable(run_tests ${TEST_SOURCE} ${MEMERNA_HEADERS})
include_directories(src)
target_link_libraries(run_tests Threads::Threads ${GTEST_LIBRARIES} memerna_lib)
