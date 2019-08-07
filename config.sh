#!/bin/sh


rm -f CMakeCache.txt

cmake \
  -C $HOME/Codes/amanzi/install/tpls/share/cmake/amanzi-tpl-config.cmake \
  ./

make VERBOSE=1
