# Download & install google test
include_guard(GLOBAL)

# Resolve the path to the git executable
find_package(Git QUIET)
if(NOT GIT_EXECUTABLE)
    message(FATAL_ERROR "error: could not resolve the path to the git executable")
endif()

# Download and update google test
set(gtest_src_directory "${PROJECT_SOURCE_DIR}/lib/google-test/") # where the library will be downloaded
if(NOT EXISTS "${gtest_src_directory}/CMakeLists.txt")
    file(REMOVE_RECURSE "${gtest_src_directory}") # just in case the folder already exists
    execute_process(COMMAND ${GIT_EXECUTABLE} clone https://github.com/google/googletest.git ${gtest_src_directory}
            RESULT_VARIABLE result)
    if(result)
        message(FATAL_ERROR "Cannot clone the remote project for googletest: ${result}")
    endif()
else() # update
    execute_process(COMMAND ${GIT_EXECUTABLE} status -s
            WORKING_DIRECTORY ${gtest_src_directory}
            OUTPUT_VARIABLE output
            ERROR_VARIABLE output
    )
    if(NOT output STREQUAL "")
       message(WARNING "The git repository `${gtest_src_directory}' contains changes and it will be not updated")
    else()
        execute_process(COMMAND ${GIT_EXECUTABLE} pull WORKING_DIRECTORY "${gtest_src_directory}" OUTPUT_QUIET)
        execute_process(COMMAND ${GIT_EXECUTABLE} checkout master WORKING_DIRECTORY "${gtest_src_directory}" OUTPUT_VARIABLE output ERROR_VARIABLE output)
        if(NOT output MATCHES "Already on")
            message(WARNING "checkout on `${gtest_src_directory}': ${output}")
        endif()
    endif()
endif()

# Prevent overriding the parent project's compiler/linker settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add the google-test project
add_subdirectory("${gtest_src_directory}" EXCLUDE_FROM_ALL)