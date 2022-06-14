include(CheckCXXCompilerFlag)

# Create a new library
add_library(common EXCLUDE_FROM_ALL
    backtrace.cpp
    cpu_topology.cpp
    database.cpp
    error.cpp
    filesystem.cpp
    math.cpp
    permutation.cpp
    profiler.cpp
    quantity.cpp
    system_compiler.cpp
    system_concurrency.cpp
    system_introspection.cpp
    timer.cpp
)

# Headers
get_target_property(path_interface interface INTERFACE_INCLUDE_DIRECTORIES)

target_link_libraries(common PUBLIC interface) # only to add transitivity of the dependency
target_include_directories(common PRIVATE "${path_interface}/lib/common") # for the internal impl~

# PAPI, for profiling
find_package(PAPI REQUIRED)
target_include_directories(common PRIVATE ${PAPI_INCLUDE_DIRS})
target_link_libraries(common PRIVATE ${PAPI_LIBRARIES})

# Dependencies
target_link_libraries(common PRIVATE libsqlite3  libbacktrace pthread dl)
target_compile_options(common PRIVATE -Wall -Werror)

# Debug flags
if((CMAKE_BUILD_TYPE STREQUAL "Debug") OR (NOT CMAKE_BUILD_TYPE))
    CHECK_CXX_COMPILER_FLAG("-fno-limit-debug-info" no_limit_debug_info)
    if(no_limit_debug_info)
        target_compile_options(common PRIVATE "-fno-limit-debug-info")
    endif()
    CHECK_CXX_COMPILER_FLAG("-fno-omit-frame-pointer" no_omit_frame_pointer)
    if(no_omit_frame_pointer)
        target_compile_options(common PRIVATE "-fno-omit-frame-pointer")
    endif()
endif()

# Numa support
find_package(Numa)
if(NUMA_FOUND)
    target_link_libraries(common PRIVATE ${NUMA_LIBRARY})
    set_property(SOURCE system_concurrency.cpp APPEND PROPERTY COMPILE_DEFINITIONS HAVE_LIBNUMA)
    # do not add -isystem /usr/include : https://gcc.gnu.org/bugzilla/show_bug.cgi?id=70129
    if(NOT ${NUMA_INCLUDE_DIR} MATCHES "/usr/include")
        set_property(SOURCE system_concurrency.cpp APPEND PROPERTY INCLUDE_DIRECTORIES ${NUMA_INCLUDE_DIR})
    endif()
endif()

