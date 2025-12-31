cmake_minimum_required(VERSION 3.14)
project(atroplex VERSION 0.0.1)
set(CMAKE_CXX_STANDARD 20)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/config.hpp.in ${CMAKE_CURRENT_SOURCE_DIR}/include/config.hpp)

# Fetch dependencies
include(FetchContent)

# cxxopts for command-line parsing
FetchContent_Declare(
        cxxopts
        GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
        GIT_TAG v3.2.1
)
FetchContent_MakeAvailable(cxxopts)

# genogrove for genomic interval data structures
FetchContent_Declare(
        genogrove
        GIT_REPOSITORY https://github.com/genogrove/genogrove.git
        GIT_TAG v0.9.0
)
FetchContent_MakeAvailable(genogrove)

# htslib for reading FASTQ/BAM files
find_package(ZLIB REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)

file(GLOB SOURCES "src/*.cpp")
add_executable(atroplex ${SOURCES})
target_include_directories(atroplex PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        ${HTSLIB_INCLUDE_DIRS}
)
target_link_libraries(atroplex
        cxxopts
        genogrove
        ${HTSLIB_LIBRARIES}
        ZLIB::ZLIB
)

