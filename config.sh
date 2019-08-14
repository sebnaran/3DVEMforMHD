#!/bin/sh


rm -rf CMakeCache.txt /build
#! Before running the cmakelists.txt we will populate the cmakecache to define a series of variables
#! Regarding the paths to the installation directories of the TPLs.
cmake \
-C $HOME/Codes/amanzi/install/tpls/share/cmake/amanzi-tpl-config.cmake \
-D CMAKE_BUILD_TYPE="Release" \
./
#!-D ENABLE_Structured:BOOL=false \
#!./

make #!VERBOSE=1
