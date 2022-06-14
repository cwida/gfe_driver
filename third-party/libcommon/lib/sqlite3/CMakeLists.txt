# Latest bundle included on 25/07/2018, v3.24.0
add_library(libsqlite3 EXCLUDE_FROM_ALL sqlite3.c sqlite3.h)
target_link_libraries(libsqlite3 dl)
if(NOT CMAKE_BUILD_TYPE) # as Debug
    target_compile_options(libsqlite3 PRIVATE -g)
endif()
set_target_properties(libsqlite3 PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR})