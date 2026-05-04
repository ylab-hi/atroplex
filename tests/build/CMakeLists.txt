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

# Hub PSI + entropy integration test: builds a grove with a known hub
# topology and verifies the .psi and .entropy columns in splicing_hubs.tsv.
add_executable(hub_psi_entropy_tests hub_psi_entropy_test.cpp)
target_link_libraries(hub_psi_entropy_tests PRIVATE
    atroplex_core
    gtest gtest_main
)

gtest_discover_tests(hub_psi_entropy_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Conserved-threshold tests: is_conserved >= semantics,
# analysis_report::set_conserved_fraction validation, and end-to-end
# verification that the conserved_segments.tsv content tracks the
# configured fraction.
add_executable(conserved_threshold_tests conserved_threshold_test.cpp)
target_link_libraries(conserved_threshold_tests PRIVATE
    atroplex_core
    gtest gtest_main
)

gtest_discover_tests(conserved_threshold_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Compact subcall tests: end-to-end disk roundtrip through subcall::compact
# (build → save → compact → reload), plus validate() failure modes
# (missing .qtx, --no-qtx opt-out, refusal to overwrite the input dir).
add_executable(compact_tests compact_test.cpp)
target_link_libraries(compact_tests PRIVATE
    atroplex_core
    gtest gtest_main
)

gtest_discover_tests(compact_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
