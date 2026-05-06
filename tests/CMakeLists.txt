include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2
)

FetchContent_MakeAvailable(googletest)
include(GoogleTest)

# Test binaries link against the shared atroplex_core static library rather
# than re-compiling the full source glob. Without this, every test binary
# triggers another N-way parallel compile of the library, which OOMs GCC on
# the CI runner.

add_subdirectory(build)
add_subdirectory(compact)
add_subdirectory(discover)
add_subdirectory(export)
add_subdirectory(inspect)
add_subdirectory(query)