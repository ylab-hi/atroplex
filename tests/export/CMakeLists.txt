# Export --conserved-fraction tests: drives subcall::export_gtf end-to-end
# with --conserved-only at various fractions, plus the (0,1] range check
# in validate().
add_executable(export_conserved_fraction_tests export_conserved_fraction_test.cpp)
target_link_libraries(export_conserved_fraction_tests PRIVATE
    atroplex_core
    gtest gtest_main
)

gtest_discover_tests(export_conserved_fraction_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Export round-trip + schema tests (#61): build a one-sample grove,
# run subcall::export_gtf end-to-end, then verify the exported file's
# schema (gene/transcript/exon lines + required attributes), exon
# coordinate fidelity, and that re-parsing the exported GTF through
# build_gff::build produces the same segment count as the original.
add_executable(export_reparse_tests export_reparse_test.cpp)
target_link_libraries(export_reparse_tests PRIVATE
    atroplex_core
    atroplex_test_utils
    gtest gtest_main
)

gtest_discover_tests(export_reparse_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
