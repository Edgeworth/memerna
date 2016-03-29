cmake_minimum_required(VERSION 3.2)

project(memerna)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FindUBSan)
include(FindASan)
include(FindMSan)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES src/main.cpp src/energy.cpp src/energy.h src/base.cpp src/base.h src/globals.cpp src/globals.h src/parsing.cpp src/parsing.h src/common.h)
add_executable(memerna ${SOURCE_FILES})
