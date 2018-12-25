# ========================================
# Author: shivanandvp 
# Email : shivanandvp.oss@gmail.com
# ========================================

# Set the minimum required cmake version
cmake_minimum_required(VERSION 2.8.9)

# Set the project name
project(NLMP_BVP)

# Set the include directories
include_directories(
***REMOVED*** include
***REMOVED*** ~/eigen
***REMOVED*** ~/eigen/unsupported
***REMOVED*** ~/eigen/unsupported/test
***REMOVED*** ~/eigen/unsupported/test/mpreal
***REMOVED*** )

# Set the source files
file(GLOB SOURCES
***REMOVED*** "src/*"
***REMOVED*** )

# Set the output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
***REMOVED*** ../bin
***REMOVED*** )

# Add an executable file to run, and specify the source files to use, to generate it
add_executable(
***REMOVED*** nlmp_bvp_test ${SOURCES}
***REMOVED*** )

# Set the directories to look for, for external libraries
LINK_DIRECTORIES(
***REMOVED*** /usr/lib
***REMOVED*** )

# Specify the external libraries to use, in order
TARGET_LINK_LIBRARIES(nlmp_bvp_test
***REMOVED*** libgmpxx.so
***REMOVED*** libgmp.so
***REMOVED*** libmpfr.so
 )
