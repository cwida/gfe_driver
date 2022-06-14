cmake_minimum_required(VERSION 3.10)

# Be verbose
set(CMAKE_VERBOSE_MAKEFILE ON)

# Prefer Clang over GCC
set(CMAKE_C_COMPILER_NAMES clang gcc icc cc)
set(CMAKE_CXX_COMPILER_NAMES clang++ g++ icpc c++ cxx)

project(common)
set(CMAKE_CXX_STANDARD 17)

# Build tests?
option(BUILD_TEST "Build the unit tests" ON)

# Search for project modules in /build-aux/
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/build-aux/")

# Helpers
if(BUILD_TEST)
    include(GoogleTest)
endif()
include(LibBacktrace)
include(SQLite3)

# Interface
add_library(interface INTERFACE)
target_include_directories(interface INTERFACE include)

# Sources, it defines the library common
add_subdirectory(src obj)

add_library(libcommon INTERFACE)
target_link_libraries(libcommon INTERFACE common)

# Tests
if(BUILD_TEST)
    add_subdirectory(test EXCLUDE_FROM_ALL)
endif()
