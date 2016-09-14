cmake_minimum_required(VERSION 3.2)

project(memerna)
set(BUILD_SHARED_LIBS OFF)

set(SAVED_BUILD_TYPE ${CMAKE_BUILD_TYPE})
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O3")
# Always have librnary compiled as release.
set(CMAKE_BUILD_TYPE Release)
add_subdirectory("extern/rnark/librnary")
# Disable warning messages for subprojects.
set_target_properties(librnary PROPERTIES COMPILE_FLAGS "-w")
set(CMAKE_BUILD_TYPE ${SAVED_BUILD_TYPE})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FindUBSan)
include(FindASan)
include(FindMSan)
include(FindThreads)

string(TOLOWER "${CMAKE_BUILD_TYPE}" SAVED_BUILD_TYPE)
if("${SAVED_BUILD_TYPE}" STREQUAL "asan" AND NOT HAVE_ADDRESS_SANITIZER)
  message(FATAL_ERROR "requested asan build but not available on current compiler")
endif()
if("${SAVED_BUILD_TYPE}" STREQUAL "msan" AND NOT HAVE_MEMORY_SANITIZER)
  message(FATAL_ERROR "requested msan build but not available on current compiler")
endif()
if(S"${AVED_BUILD_TYPE}" STREQUAL "ubsan" AND NOT HAVE_UNDEFINED_BEHAVIOR_SANITIZER)
  message(FATAL_ERROR "requested ubsan build but not available on current compiler")
endif()

set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -Wshadow \
  -Wconversion -Wfloat-equal -Wformat-security -Wno-sign-conversion -pedantic")
set(CMAKE_CXX_FLAGS_ASAN "${CMAKE_CXX_FLAGS_ASAN} -march=native -O3")
set(CMAKE_CXX_FLAGS_MSAN "${CMAKE_CXX_FLAGS_MSAN}  -march=native -O3")
set(CMAKE_CXX_FLAGS_UBSAN "${CMAKE_CXX_FLAGS_UBSAN} -march=native -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -march=native -O3")
file(GLOB MEMERNA_SOURCE
    "src/*.cpp" "src/*.h" "src/energy/*.cpp" "src/energy/*.h" "src/fold/*.cpp"
    "src/fold/*.h")
file(GLOB BRIDGE_SOURCE "src/bridge/*.cpp" "src/bridge/*.h")
file(GLOB TEST_SOURCE "tests/*.cpp" "tests/*.h")

include_directories(src)
# Suppress warnings.
include_directories(SYSTEM extern/rnark/librnary/include)
include_directories(SYSTEM extern/rnark/librnary/lib)

add_library(memerna ${MEMERNA_SOURCE})
add_library(bridge ${BRIDGE_SOURCE})

add_executable(efn src/programs/efn.cpp)
add_executable(fold src/programs/fold.cpp)
add_executable(fuzz src/programs/fuzz.cpp)
add_executable(subopt src/programs/subopt.cpp)
add_executable(harness src/programs/harness.cpp)

add_subdirectory(extern/gtest)
enable_testing()
add_executable(run_tests ${TEST_SOURCE})

target_link_libraries(run_tests memerna Threads::Threads gtest)
target_link_libraries(efn memerna)
target_link_libraries(fold memerna)
target_link_libraries(fuzz bridge memerna librnary)
target_link_libraries(subopt memerna)
target_link_libraries(harness bridge memerna librnary)
