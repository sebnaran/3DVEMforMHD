#!/bin/sh


rm -f CMakeCache.txt
#! Before running the cmakelists.txt we will populate the cmakecache to define a series of variables
#! Regarding the paths to the installation directories of the TPLs.
cmake \
-C $HOME/Codes/amanzi/install/tpls/share/cmake/amanzi-tpl-config.cmake \
./

make VERBOSE=1