# Compact subcall tests: end-to-end disk roundtrip through subcall::compact
# (build → save → compact → reload), plus validate() failure modes
# (missing .qtx, --no-qtx opt-out, refusal to overwrite the input dir).
add_executable(compact_tests compact_test.cpp)
target_link_libraries(compact_tests PRIVATE
    atroplex_core
    gtest gtest_main
)

gtest_discover_tests(compact_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})