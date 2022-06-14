# Download and build libbacktrace. It creates the target `libbacktrace' to be used.
include_guard(GLOBAL)
include(CheckCXXSourceCompiles)
include(ExternalProject)

# Dependency on cxxabi.h
#Find_Path(CXXABI_INCLUDE_DIR cxxabi.h) # This tends to fail
check_cxx_source_compiles("
#include <cxxabi.h>
int main(int argc, char* argv[]){
    char * type; int status;
    char * r = abi::__cxa_demangle(type, 0, 0, &status);
    return 0;
}\n" CXXABI)
if(NOT CXXABI)
    message(FATAL_ERROR "Missing mandatory dependency on cxxabi.h")
endif()


# Download & build libbacktrace
ExternalProject_Add(libbacktrace_external
        GIT_REPOSITORY "https://github.com/ianlancetaylor/libbacktrace" # upstream repository
        GIT_TAG "master"
        SOURCE_DIR "${PROJECT_SOURCE_DIR}/lib/backtrace" # where to download the library
        BINARY_DIR "${CMAKE_BINARY_DIR}/lib/backtrace/build" # where the library is going to be built
        PREFIX "${CMAKE_BINARY_DIR}/lib/backtrace/cmake" # cmake logs and other rubbish
        INSTALL_DIR "${CMAKE_BINARY_DIR}/lib/backtrace" # not sure?
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --disable-shared --enable-static
        BUILD_COMMAND make -j # the command to build the library
        BUILD_BYPRODUCTS "${CMAKE_BINARY_DIR}/lib/backtrace/lib/libbacktrace.a"
        UPDATE_COMMAND "" # hack, avoid to rebuild the whole library once built
        PATCH_COMMAND ""
)

# Create the dependency
add_library(libbacktrace STATIC IMPORTED)
add_dependencies(libbacktrace libbacktrace_external)
# ensure that the include directory exists even before the libbacktrace is created
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/lib/backtrace/include")
set_target_properties(libbacktrace PROPERTIES
        IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/lib/backtrace/lib/libbacktrace.a"
        INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/lib/backtrace/include")
