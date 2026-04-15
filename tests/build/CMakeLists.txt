# Absorption rule tests
add_executable(absorption_tests absorption_test.cpp)
target_link_libraries(absorption_tests PRIVATE
    atroplex_core
    gtest gtest_main
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
add_executable(serialization_tests serialization_test.cpp)
target_link_libraries(serialization_tests PRIVATE
    atroplex_core
    gtest gtest_main
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

# Builder pipeline tests: builder::build_from_samples, remove_tombstones,
# build_counters aggregation, build_summary::write_summary
add_executable(builder_pipeline_tests builder_pipeline_test.cpp)
target_link_libraries(builder_pipeline_tests PRIVATE
    atroplex_core
    gtest gtest_main
)

gtest_discover_tests(builder_pipeline_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
