add_executable(absorption_tests absorption_test.cpp ${ATROPLEX_LIB_SOURCES})
target_include_directories(absorption_tests PRIVATE
    "${CMAKE_SOURCE_DIR}/include"
    ${HTSLIB_INCLUDE_DIRS}
)
target_link_libraries(absorption_tests
    gtest gtest_main
    genogrove
    cxxopts
    ${HTSLIB_LIBRARIES}
    ZLIB::ZLIB
)

# Copy fixtures at build time (not just configure time) so new files are always picked up
add_custom_command(TARGET absorption_tests POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/fixtures
        ${CMAKE_CURRENT_BINARY_DIR}/fixtures
)
target_compile_definitions(absorption_tests PRIVATE
    TEST_FIXTURE_DIR="${CMAKE_CURRENT_BINARY_DIR}/fixtures"
)

gtest_discover_tests(absorption_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Serialization roundtrip tests
add_executable(serialization_tests serialization_test.cpp ${ATROPLEX_LIB_SOURCES})
target_include_directories(serialization_tests PRIVATE
    "${CMAKE_SOURCE_DIR}/include"
    ${HTSLIB_INCLUDE_DIRS}
)
target_link_libraries(serialization_tests
    gtest gtest_main
    genogrove
    cxxopts
    ${HTSLIB_LIBRARIES}
    ZLIB::ZLIB
)

add_custom_command(TARGET serialization_tests POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/fixtures
        ${CMAKE_CURRENT_BINARY_DIR}/fixtures
)
target_compile_definitions(serialization_tests PRIVATE
    TEST_FIXTURE_DIR="${CMAKE_CURRENT_BINARY_DIR}/fixtures"
)

gtest_discover_tests(serialization_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})