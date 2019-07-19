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

# Determine the build type. If none is specified, assume release version.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

set(CMAKE_CXX_FLAGS "-Wall -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# This is a thing that will help with emacs code completion
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# This was the default boost on my Ubuntu 18 machine, I'm not tied to it in any way -- JSF
# To add other boost libraries, just keep adding to the list following COMPONENTS
find_package(Boost 1.65.0 COMPONENTS program_options REQUIRED)

# NLOPT setup. Download from: https://github.com/stevengj/nlopt/archive/v2.6.1.tar.gz
# Extract and follow directions here: https://nlopt.readthedocs.io/en/latest/NLopt_Installation/
find_package(NLopt COMPONENTS program_options REQUIRED)

# Eigen setup. Download from: http://eigen.tuxfamily.org/index.php?title=Main_Page#Download
# Extract and follow INSTALL directions for global install
# Make sure the EIGEN3_INCLUDE_DIR has been set to the directory containing the 'Eigen' folder in the download:
# "export EIGEN3_INCLUDE_DIR=/path/to/eigenfolder/" is nice to put in your VM's ~/.bashrc
set( EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
if( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
else()
    MESSAGE(STATUS "Using Eigen directory: ${EIGEN3_INCLUDE_DIR}")
endif()

set(EIGEN3_INC_DIR $ENV{EIGEN3_PATH})

# Set any specific compile options you need here:
add_compile_options(-std=c++11 -Wall)

# Set up the libgtoc library that executables will link with
# If you add a source file (.cpp/.cc), add it to the
set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(TEST_DIR ${CMAKE_SOURCE_DIR}/test)
set(BIN_DIR ${CMAKE_SOURCE_DIR}/bin)
set(INC_DIR ${CMAKE_SOURCE_DIR}/inc)
include_directories(${INC_DIR} ${Boost_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${NLOPT_INCLUDE_DIR})
set(SOURCES ${SOURCE_DIR}/propagate.cpp ${SOURCE_DIR}/star.cpp ${SOURCE_DIR}/events.cpp ${SOURCE_DIR}/transfer.cpp)
add_library(gtoc ${SOURCES})
set(external_link_libraries ${Boost_LIBRARIES} ${NLOPT_LIBRARIES})

set(EXECUTABLE_OUTPUT_PATH ${BIN_DIR})

# Do the following for each executable
add_executable(transfer_test ${CMAKE_SOURCE_DIR}/test/transfer_test.cpp)
target_link_libraries(transfer_test ${external_link_libraries} gtoc nlopt)

add_executable(num_jac_test ${CMAKE_SOURCE_DIR}/test/num_jac_test.cpp)
target_link_libraries(num_jac_test ${external_link_libraries} gtoc nlopt)

add_executable(star_test ${TEST_DIR}/star_test.cpp)
target_link_libraries(star_test ${external_link_libraries} gtoc)

add_executable(event_test ${TEST_DIR}/event_test.cpp)
target_link_libraries(event_test ${external_link_libraries} gtoc)

add_executable(prop_test ${TEST_DIR}/prop_test.cpp)
target_link_libraries(prop_test ${external_link_libraries} gtoc)

add_executable(merit_function_test ${TEST_DIR}/merit_function_test.cpp)
target_link_libraries(merit_function_test ${external_link_libraries} gtoc)
