project(SourcesFiles C CXX)

include_directories(${SourceFiles_SOURCE_DIRS})
include_directories(../../amanzi/install/amanzi/include)
include_directories(../../amanzi/install/tpls/include)
include_directories(../..//amanzi/install/tpls/include/UnitTest++)
include_directories(../../amanzi/install/tpls/trilinos-12-12-1/include)
include_directories(../../amanzi/install/tpls/petsc-3.8.2/include/petsc/mpiuni)
include_directories(../../amanzi/install/tpls/petsc-3.8.2/include)

add_executable(MHD main.cc)
