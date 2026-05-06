# Export --conserved-fraction tests: drives subcall::export_gtf end-to-end
# with --conserved-only at various fractions, plus the (0,1] range check
# in validate().
add_executable(export_conserved_fraction_tests export_conserved_fraction_test.cpp)
target_link_libraries(export_conserved_fraction_tests PRIVATE
    atroplex_core
    gtest gtest_main
)

gtest_discover_tests(export_conserved_fraction_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
