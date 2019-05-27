# GTOCX

## Building the project
Currently, things are collected, compiled and linked using CMake.
All you need to do is build a "build" directory, cd into it, and call
cmake on its parent directory. When cmake is finished, run make and 
executables will compile to the bin directory. To add new executables,
follow the examples on the bottom of the CMakeLists.txt file.

build.sh will build the project, or rerun cmake and make
clean.sh will reset the project build directory to github state.
