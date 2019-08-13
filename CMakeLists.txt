cmake_minimum_required(VERSION 3.10)
#With this command I am enabling the variables VEM_SOURCE_DIR which specifies the path to
#headers and footers and the VEM_BINARY_DIR which is where all binary executables will reside.
#These are defaulted to ./ unless otherwise specified.
project(VEM LANGUAGES CXX)

#Here we will specify that we want all binary files to be placed in ./build/bin
#set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)

#Now we specify that we want all executables to be placed in the binary directory defined above.
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR})
#set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

#We will need the addresses of several variables throughout amanzi, so we begin by speficying the 
#module path
#set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../amanzi/build/tpls/trilinos/trilinos-12-12-1-build/packages/teuchos)

#This module will populate the variables with the list of dependecies of the Teuchos libraries.
#include(TeuchosConfig)

add_subdirectory(src)
