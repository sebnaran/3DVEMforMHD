cmake_minimum_required(VERSION 3.6)
#With this command I am enabling the variables VEM_SOURCE_DIR which specifies the path to
#headers and footers and the VEM_BINARY_DIR which is where all binary executables will reside.
#These are defaulted to ./ unless otherwise specified.
project(VEM LANGUAGES CXX)
#Here we will specify that we want all binary files to be placed in ./build/bin
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)
#Now we specify that we want all executables to be placed in the binary directory defined above.
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

add_subdirectory(src)