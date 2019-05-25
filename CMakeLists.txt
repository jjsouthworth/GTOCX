# Author(s): jsf
# cmake for GTOC X project.
# If you don't have cmake: 'sudo apt install cmake'

# Building instructions:
# cd $GTOC_TOP_DIRECTORY; mkdir build && cd build; cmake ..; make


#First few lines are setup, some project specifications, etc
cmake_minimum_required (VERSION 3.10)
project(GTOCX)
set(GTOCX_VERSION_MAJOR_0)
set(GTOCX_VERSION_MINOR_1)

# This was the default boost on my Ubuntu 18 machine, I'm not tied to it in any way -- JSF
# To add other boost libraries, just keep adding to the list following COMPONENTS
find_package(Boost 1.65.0 COMPONENTS program_options REQUIRED)

# Set up the libgtoc library that executables will link with
# If you add a source file (.cpp/.cc), add it to the 
include_directories(inc ${Boost_INCLUDE_DIR})
set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(SOURCES ${SOURCE_DIR}/propagate.cpp ${SOURCE_DIR}/star.cpp ${SOURCE_DIR}/vessels.cpp)
add_library(gtoc ${SOURCES})

# Do the following for each executable
add_executable(star_test ${CMAKE_SOURCE_DIR}/test/star_test.cpp)
target_link_libraries(star_test ${Boost_LIBRARIES} gtoc)
