project(SourcesFiles C CXX)

include_directories(${SourceFiles_SOURCE_DIRS})
include_directories(../../amanzi/install/amanzi/include)
include_directories(../../amanzi/install/tpls/include)
include_directories(../..//amanzi/install/tpls/include/UnitTest++)
include_directories(../../amanzi/install/tpls/trilinos-12-12-1/include)
include_directories(../../amanzi/install/tpls/petsc-3.8.2/include/petsc/mpiuni)
include_directories(../../amanzi/install/tpls/petsc-3.8.2/include)
include_directories(../../amanzi/src/operators/test)
include_directories(${LD_LIBRARY_PATH}/../include)

#The definition of the addresses found in the "target_link_libraries" command are 
#found in $HOME/Codes/amanzi/install/tpls/share/cmake/amanzi-tpl-config.cmake


link_directories(../../amanzi/install/amanzi/lib 
                 ../../amanzi/install/tpls/lib
                 ../../amanzi/install/tpls/trilinos-12-12-1/lib
                 ../../amanzi/install/tpls/SEACAS/lib
                 /usr/lib/x86_64-linux-gnu/) 

#########################################################################################################
# #Attempt number 1
 set(Amanzi_libs mesh mesh_factory mstk_mesh simple_mesh geometry 
                 operators solvers data_structures whetstone output atk error_handling)

set(Teuchos_libs teuchoskokkoscomm teuchoskokkoscompat teuchosremainder teuchosnumerics teuchoscomm
                 teuchosparameterlist teuchoscore kokkoscore)

 set(Epetra_libs epetra epetraext teuchoskokkoscomm teuchoskokkoscompat teuchosremainder teuchosnumerics 
                 teuchoscomm teuchosparameterlist teuchoscore kokkoscore)

set(Epetra_interface_libs teuchoskokkoscomm teuchoskokkoscompat teuchosremainder teuchosnumerics
                          teuchoscomm teuchosparameterlist teuchoscore liblapack.so liblapack.so liblapack.so)

set(Nox_libs locaepetra localapack loca noxepetra noxlapack nox belostpetra belosepetra belos ml ifpack
             amesos aztecoo epetraext trilinosss tpetraext tpetrainout tpetra tpetraclassiclinalg
             tpetraclassicnodeapi tpetraclassic triutils zoltan epetra kokkoskernels teuchoskokkoscomm
             teuchoskokkoscompat teuchosremainder teuchosnumerics teuchoscomm teuchosparameterlist
             teuchoscore kokkosalgorithms kokkoscontainers kokkoscore)

set(Nox_interface_libs belostpetra belosepetra belos ml ifpack amesos aztecoo epetraext epetra
                        teuchoskokkoscomm teuchoskokkoscompat teuchosremainder teuchosnumerics teuchoscomm
                        teuchosparameterlist teuchoscore liblapack.so libblas.so libblas.so)

set (Boost_libs boost_system boost_filesystem)

set(HDF5_libs hdf5_hl hdf5 hdf5_hl libz.so m dl)

set(TPLS_meshes mstk metis exodus netcdf)

# # set(amanzi_libs atk chemistry_pk data_structures dbg dummy energy eos error_handling flow functions 
# #                 geochembase geochemrxns geochemsolvers geochemutil geometry interface interface_s mesh
# #                 mesh_audit mesh_factory mesh_functions mesh_logical mpc mpc_tree mstk_mesh navier_stokes
# #                 operators output pks simple_mesh solvers state structured_amanzi time_integration transport
# #                 udriv whetstone)

#  set(trilinos_libs amesos amesos2 aztecoo belos belosepetra belostpetra epetra epetraext gtest ifpack 
#                    kokkosalgorithms kokkoscontainers kokkoscore kokkoskernels loca locaepetra localapack
#                    ml nox noxepetra noxlapack teuchoscomm teuchoscore teuchoskokkoscomm teuchoskokkoscompat
#                    teuchosnumerics teuchosparameterlist teuchosremainder tpetra tpetraclassic tpetraclassiclinalg
#                    tpetraclassicnodeapi tpetraext tpetrainout trilinosss triutils zoltan)

# This builds the Test executable.
add_executable(TestMHD TestMain.cc)

target_link_libraries(TestMHD ${Amanzi_libs} ${TPLS_meshes} ${Epetra_libs} ${Epetra_interface_libs} ${Teuchos_libs}  
                              ${Nox_libs} ${Nox_interface_libs} HYPRE superlu_dist superlu parmetis metis
                              UnitTest++ epetra ${Boost_libs} ${HDF5_libs})

#target_link_libraries(TestMHD ${Amanzi_libs} ${Epetra_libs} UnitTest++ ${HDF5_libs})
#########################################################################################################



#Attempt number 2



#add_executable(TestMHD TestMain.cc)

#target_link_libraries(TestMHD UnitTest++)









































# target_link_libraries(TestMHD ${Amanzi_libs} ${Epetra_libs} ${Epetra_interface_libs} ${Teuchos_libs}  
#                               ${Nox_libs} ${Nox_interface_libs} UnitTest++ ${HDF5_libs})

#target_link_libraries(TestMHD ${UnitTest_DIR}/lib/libUnitTest++.a)
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
