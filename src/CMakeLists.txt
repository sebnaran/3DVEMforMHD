project(SourcesFiles C CXX)

include_directories(${SourceFiles_SOURCE_DIRS})
include_directories(../../amanzi/install/amanzi/include)
include_directories(../../amanzi/install/tpls/include)
include_directories(../..//amanzi/install/tpls/include/UnitTest++)
include_directories(../../amanzi/install/tpls/trilinos-12-12-1/include)
include_directories(../../amanzi/install/tpls/petsc-3.8.2/include/petsc/mpiuni)
include_directories(../../amanzi/install/tpls/petsc-3.8.2/include)
include_directories(../../amanzi/src/operators/test)
#The definition of the addresses found in the "target_link_libraries" command are 
#found in $HOME/Codes/amanzi/install/tpls/share/cmake/amanzi-tpl-config.cmake

#This builds the Test executable.
add_executable(TestMHD TestMain.cc)

# target_link_libraries(TestMHD PUBLIC /home/sebastian/Codes/amanzi/install/amanzi/lib/liboperators.a
#                               PUBLIC ${UnitTest_DIR}/lib/libUnitTest++.a 
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libteuchoscore.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libteuchosnumerics.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libteuchosremainder.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libteuchosparameterlist.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libteuchoskokkoscompat.a 
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libteuchoskokkoscomm.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libamesos2.a 
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libamesos.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libaztecoo.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libbelos.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libbelosepetra.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libbelostpetra.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libifpack.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libepetra.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libepetraext.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libgtest.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libkokkosalgorithms.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libkokkoscontainers.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libkokkoscore.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libkokkoskernels.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libloca.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/liblocaepetra.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/liblocalapack.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libml.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libnox.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libnoxepetra.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libnoxlapack.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libteuchoscomm.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libtpetra.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libtpetraclassic.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libtpetraclassiclinalg.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libtpetraclassicnodeapi.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libtpetraext.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libtpetrainout.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libtrilinosss.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libtriutils.a
#                               PUBLIC ${Trilinos_INSTALL_PREFIX}/lib/libzoltan.a
#                               PUBLIC /home/sebastian/Codes/amanzi/install/amanzi/lib/libmesh_factory.a
#                               PUBLIC /home/sebastian/Codes/amanzi/install/amanzi/lib/libmesh.a
#                               PUBLIC /home/sebastian/Codes/amanzi/install/amanzi/lib/liboperators.a
#                               PUBLIC /home/sebastian/Codes/amanzi/install/amanzi/lib/libwhetstone.a
#                               /usr/lib/x86_64-linux-gnu/libdl.so)



#This builds the main code
#add_executable(MHD main.cc)
#target_link_libraries(MHD ${Trilinos_INSTALL_PREFIX}/lib/libteuchoscore.a)