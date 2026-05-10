add_executable(query_tests query_classification_test.cpp)
target_link_libraries(query_tests PRIVATE
    atroplex_core
    gtest gtest_main
)

add_custom_command(TARGET query_tests POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/fixtures
        ${CMAKE_CURRENT_BINARY_DIR}/fixtures
)
target_compile_definitions(query_tests PRIVATE
    TEST_FIXTURE_DIR="${CMAKE_CURRENT_BINARY_DIR}/fixtures"
)

gtest_discover_tests(query_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Closest-match (flanking-based nearest neighbor) tests — fixture-free,
# builds the grove inline.
add_executable(closest_match_tests closest_match_test.cpp)
target_link_libraries(closest_match_tests PRIVATE
    atroplex_core
    gtest gtest_main
)
gtest_discover_tests(closest_match_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Writer-level tests for query::write_classification — drives the writer
# directly with a synthetic vector<query_result>, asserts column schema
# with/without --nearest. Fixture-free; reads the temp output file back.
add_executable(write_classification_tests write_classification_test.cpp)
target_link_libraries(write_classification_tests PRIVATE
    atroplex_core
    gtest gtest_main
)
gtest_discover_tests(write_classification_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
