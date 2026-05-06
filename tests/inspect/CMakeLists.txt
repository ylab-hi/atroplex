# Hub PSI + entropy integration tests: builds a grove with a known hub
# topology and verifies splicing_hubs.tsv (PSI / entropy / branches columns,
# default-10 and configurable-threshold boundary cases) plus the
# branch_details.tsv schema and content.
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
