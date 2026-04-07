include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2
)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

# Collect atroplex source files (excluding main.cpp)
file(GLOB ATROPLEX_LIB_SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp" "${CMAKE_SOURCE_DIR}/src/subcall/*.cpp")
list(FILTER ATROPLEX_LIB_SOURCES EXCLUDE REGEX "main\\.cpp$")

add_subdirectory(build)