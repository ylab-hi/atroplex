add_executable(discover_tests sqanti_category_test.cpp ${ATROPLEX_LIB_SOURCES})
target_include_directories(discover_tests PRIVATE
    "${CMAKE_SOURCE_DIR}/include"
    ${HTSLIB_INCLUDE_DIRS}
)
target_link_libraries(discover_tests
    gtest gtest_main
    genogrove
    cxxopts
    ${HTSLIB_LIBRARIES}
    ZLIB::ZLIB
)

# Copy fixtures at build time so new files are always picked up
add_custom_command(TARGET discover_tests POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/fixtures
        ${CMAKE_CURRENT_BINARY_DIR}/fixtures
)
target_compile_definitions(discover_tests PRIVATE
    TEST_FIXTURE_DIR="${CMAKE_CURRENT_BINARY_DIR}/fixtures"
)

gtest_discover_tests(discover_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})