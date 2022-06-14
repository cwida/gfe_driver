# Based on https://github.com/HewlettPackard/foedus_code/blob/master/foedus-core/cmake/FindNuma.cmake

# Find the numa policy library.
# Output variables:
#  NUMA_INCLUDE_DIR : e.g., /usr/include/.
#  NUMA_LIBRARY     : Library path of numa library
#  NUMA_FOUND       : True if found.

find_path(NUMA_INCLUDE_DIR NAMES numa.h)
find_library(NUMA_LIBRARY NAMES numa)

# Handle the QUIETLY and REQUIRED arguments and set SQLITE3_FOUND to TRUE if all listed variables are TRUE.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Numa DEFAULT_MSG NUMA_LIBRARY NUMA_INCLUDE_DIR)
mark_as_advanced(NUMA_INCLUDE_DIR NUMA_LIBRARY)