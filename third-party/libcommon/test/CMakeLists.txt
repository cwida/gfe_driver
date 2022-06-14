file(GLOB SOURCES *.cpp)
add_executable(libcommon_tests ${SOURCES})
target_link_libraries(libcommon_tests libcommon gtest_main)
add_test(NAME libcommon_tests COMMAND libcommon_tests)