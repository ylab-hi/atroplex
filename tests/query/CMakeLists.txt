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
