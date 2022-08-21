# ========================================
# Author: shivanandvp 
# Email : shivanandvp.oss@gmail.com
# ========================================

# Set the minimum required cmake version
cmake_minimum_required(VERSION 3.13)

# Set the project name
project(
    nlmpBVP       
        VERSION 1.0.0 
        LANGUAGES CXX 
    )

# create a list of files "execSources" to include in the target executable
file(
    GLOB_RECURSE execSources 
        src/*.cpp    
    )

# create a list of files "execHeaders" to include in the target executable
file(
    GLOB_RECURSE execHeaders 
        include/*.hpp                  # add the header files
    )

# Add a target executable file to run, and specify the source files to use, to generate it
add_executable(
    test_nlmpBVP   
        ${execSources}
        ${execHeaders}
    )

target_include_directories(
    test_nlmpBVP
        PUBLIC
            include
        PRIVATE
            /usr/include/eigen3
            /usr/include/eigen3/unsupported
)

# Set the directories to look for, for external libraries
target_link_directories(
    test_nlmpBVP
        PRIVATE
            /usr/lib
    )

# Specify the external libraries to use, in order
target_link_libraries(
    test_nlmpBVP
        PRIVATE
            libgmpxx.so
            libgmp.so
            libmpfr.so
)

# Set the output directory
set_target_properties(
    test_nlmpBVP
        PROPERTIES
            CMAKE_RUNTIME_OUTPUT_DIRECTORY
                ../bin
    )
