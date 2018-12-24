# ========================================
# Author: Shivanand Pattanshetti 
# Email : shivanand.pattanshetti@gmail.com
# ========================================

# Set the minimum required cmake version
cmake_minimum_required (VERSION 2.8.9)

# Set the project name
project(NLMP_BVP)

# Set the include directories
include_directories(
    include
    ~/eigen
    ~/eigen/unsupported
    ~/eigen/unsupported/test
    ~/eigen/unsupported/test/mpreal
    )

# Set the source files
file(GLOB SOURCES
    "src/*"
    )

# Set the output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ../bin
    )

# Set the directories to look for, for external libraries
LINK_DIRECTORIES(
    /usr/lib
    )

# Add an executable file to run, and specify the source files to use, to generate it
add_executable(
    nlmp_bvp_test ${SOURCES}
    )

# Specify the external libraries to use, in order
TARGET_LINK_LIBRARIES(nlmp_bvp_test
    libgmpxx.so
    libgmp.so
    libmpfr.so
 )
