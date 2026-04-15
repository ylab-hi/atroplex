cmake_minimum_required(VERSION 3.14)
project(atroplex VERSION 0.0.1)
set(CMAKE_CXX_STANDARD 20)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/config.hpp.in ${CMAKE_CURRENT_SOURCE_DIR}/include/config.hpp)

# Save user's BUILD_TESTING preference, then disable for subdependencies
if(NOT DEFINED ATROPLEX_BUILD_TESTING)
    set(ATROPLEX_BUILD_TESTING ON)
endif()
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

# Fetch dependencies
include(FetchContent)

# cxxopts for command-line parsing
FetchContent_Declare(
        cxxopts
        GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
        GIT_TAG v3.3.1
)
FetchContent_MakeAvailable(cxxopts)

# genogrove for genomic interval data structures
FetchContent_Declare(
        genogrove
        GIT_REPOSITORY https://github.com/genogrove/genogrove.git
        GIT_TAG v0.21.0
)
FetchContent_MakeAvailable(genogrove)

# htslib for reading FASTQ/BAM files
find_package(ZLIB REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)

file(GLOB SOURCES "src/*.cpp")
file(GLOB SUBCALL_SOURCES "src/subcall/*.cpp")

# Library target containing every source except main.cpp. The atroplex
# executable and every test binary link against this static library, so
# the library sources are compiled exactly once. Without this, each test
# binary re-compiles the full source glob and we blow past the CI
# runner's RAM budget during parallel builds.
list(FILTER SOURCES EXCLUDE REGEX "main\\.cpp$")
add_library(atroplex_core STATIC ${SOURCES} ${SUBCALL_SOURCES})
target_include_directories(atroplex_core PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        ${HTSLIB_INCLUDE_DIRS}
)
target_link_libraries(atroplex_core PUBLIC
        cxxopts
        genogrove
        ${HTSLIB_LIBRARIES}
        ZLIB::ZLIB
)

add_executable(atroplex src/main.cpp)
target_link_libraries(atroplex PRIVATE atroplex_core)

# Restore testing if user didn't explicitly disable it and tests/ exists
if(ATROPLEX_BUILD_TESTING AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    set(BUILD_TESTING ON CACHE BOOL "" FORCE)
    enable_testing()
    add_subdirectory(tests)
endif()

