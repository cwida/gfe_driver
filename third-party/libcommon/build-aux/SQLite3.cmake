# Sqlite3 bundle
set(SQLITE3_SYSTEM "ON" CACHE BOOL "Use the installation of SQLite3 from the system, if present.")
mark_as_advanced(SQLITE3_SYSTEM)

if(SQLITE3_SYSTEM)
    find_package(SQLite3)
endif()

if(NOT SQLITE3_SYSTEM OR NOT SQLITE3_FOUND)
    message(STATUS "Using the bundled SQLite3 library")
    add_subdirectory("${PROJECT_SOURCE_DIR}/lib/sqlite3/")
endif()