cmake_minimum_required (VERSION 2.8.9)

project(NLMP_BVP)

include_directories(include /usr/include/eigen3)
file(GLOB SOURCES "src/*")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)

add_executable(nlmp_bvp_test ${SOURCES})
