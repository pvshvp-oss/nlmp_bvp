cmake_minimum_required (VERSION 2.8.9)

project(NLMP_BVP)

include_directories(
    include
    ~/eigen
    ~/eigen/unsupported
    ~/eigen/unsupported/test
    ~/eigen/unsupported/test/mpreal
)
file(GLOB SOURCES "src/*")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)

LINK_DIRECTORIES(/usr/lib)

add_executable(nlmp_bvp_test ${SOURCES})

TARGET_LINK_LIBRARIES(nlmp_bvp_test
    libgmpxx.so
    libgmp.so
    libmpfr.so
 )
