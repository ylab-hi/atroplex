# Absorption rule tests
add_executable(absorption_tests absorption_test.cpp)
target_link_libraries(absorption_tests PRIVATE
    atroplex_core
    atroplex_test_utils
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
    atroplex_test_utils
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
    atroplex_test_utils
    gtest gtest_main
)

gtest_discover_tests(builder_pipeline_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Utility helpers: is_main_chromosome scaffold filter edge cases.
# Pure unit tests, no GTF fixtures needed.
add_executable(utility_tests utility_test.cpp)
target_link_libraries(utility_tests PRIVATE
    atroplex_core
    gtest gtest_main
)

gtest_discover_tests(utility_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Quantification sidecar (.qtx) round-trip tests.
# Tests Writer/Reader in isolation, no build pipeline involvement.
add_executable(quant_sidecar_tests quant_sidecar_test.cpp)
target_link_libraries(quant_sidecar_tests PRIVATE
    atroplex_core
    gtest gtest_main
)

gtest_discover_tests(quant_sidecar_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Sample manifest parser edge cases (#61): required-column enforcement,
# header / field whitespace tolerance, '.' empty marker, auto-inferred
# id and replicate group, expression-attribute splitting.
add_executable(sample_manifest_tests sample_manifest_test.cpp)
target_link_libraries(sample_manifest_tests PRIVATE
    atroplex_core
    gtest gtest_main
)

gtest_discover_tests(sample_manifest_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Malformed GTF input tests (#61): pins current behavior on degenerate
# inputs (empty file, missing gene_id / transcript_id, non-existent
# path, truncated last line). Most cases silently `continue` past the
# bad record; tests are regression markers for any future warning
# counters or strict-mode rejection.
add_executable(malformed_gtf_tests malformed_gtf_test.cpp)
target_link_libraries(malformed_gtf_tests PRIVATE
    atroplex_core
    atroplex_test_utils
    gtest gtest_main
)

gtest_discover_tests(malformed_gtf_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
