# ========================================
# Author: shivanandvp 
# Email : shivanandvp.oss@gmail.com
# ========================================

# Set the minimum required cmake version
cmake_minimum_required(VERSION 3.13)

# Set the project name
project(
***REMOVED*** nlmpBVP***REMOVED******REMOVED*** 
***REMOVED******REMOVED***  VERSION 1.0.0 
***REMOVED******REMOVED***  LANGUAGES CXX 
***REMOVED*** )

# create a list of files "execSources" to include in the target executable
file(
***REMOVED*** GLOB_RECURSE execSources 
***REMOVED******REMOVED***  src/*.cpp***REMOVED*** 
***REMOVED*** )

# create a list of files "execHeaders" to include in the target executable
file(
***REMOVED*** GLOB_RECURSE execHeaders 
***REMOVED******REMOVED***  include/*.hpp***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***# add the header files
***REMOVED*** )

# Add a target executable file to run, and specify the source files to use, to generate it
add_executable(
***REMOVED*** test_nlmpBVP***REMOVED***
***REMOVED******REMOVED***  ${execSources}
***REMOVED******REMOVED***  ${execHeaders}
***REMOVED*** )

# Set the include directories***REMOVED***
target_include_directories(
***REMOVED*** test_nlmpBVP
***REMOVED******REMOVED***  PUBLIC
***REMOVED******REMOVED******REMOVED******REMOVED***include
***REMOVED******REMOVED***  PRIVATE
***REMOVED******REMOVED******REMOVED******REMOVED***~/eigen
***REMOVED******REMOVED******REMOVED******REMOVED***~/eigen/unsupported
***REMOVED******REMOVED******REMOVED******REMOVED***~/eigen/unsupported/test
***REMOVED******REMOVED******REMOVED******REMOVED***~/eigen/unsupported/test/mpreal
)

# Set the directories to look for, for external libraries
target_link_directories(
***REMOVED*** test_nlmpBVP
***REMOVED******REMOVED***  PRIVATE
***REMOVED******REMOVED******REMOVED******REMOVED***/usr/lib
***REMOVED*** )

# Specify the external libraries to use, in order
target_link_libraries(
***REMOVED*** test_nlmpBVP
***REMOVED******REMOVED***  libgmpxx.so
***REMOVED******REMOVED***  libgmp.so
***REMOVED******REMOVED***  libmpfr.so
)

# Set the output directory
set_target_properties(
***REMOVED*** test_nlmpBVP
***REMOVED******REMOVED***  PROPERTIES
***REMOVED******REMOVED******REMOVED******REMOVED***CMAKE_RUNTIME_OUTPUT_DIRECTORY
***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** ../bin
***REMOVED*** )


