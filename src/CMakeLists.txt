project(SourcesFiles C CXX)

#set(AMANZI /Users/lipnikov/ASCEM/amanzi-install)
#set(TPLS /Users/lipnikov/ASCEM/TPLs-clang-0.95.7-ompi-dso)
set(AMANZI ../../amanzi/install/amanzi)
set(TPLS ../../amanzi/install/tpls)
include_directories(${SourceFiles_SOURCE_DIRS})
include_directories(${AMANZI}/include)
include_directories(${TPLS}/include)
include_directories(${TPLS}/include/UnitTest++)
include_directories(${TPLS}/trilinos-12-14-0bf557717e-Jun17/include)
include_directories(${TPLS}/petsc-3.8.2/include)
link_directories(${AMANZI}/lib 
                 ${TPLS}/lib
                 ${TPLS}/trilinos-12-14-0bf557717e-Jun17/lib
                 ${TPLS}/SEACAS/lib
                 /usr/lib/x86_64-linux-gnu/)

#########################################################################################################
set(Amanzi_libs mesh mesh_factory mesh_extracted mstk_mesh simple_mesh geometry 
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

set(Boost_libs boost_system boost_filesystem)

set(HDF5_libs hdf5_hl hdf5 hdf5_hl libz.so m dl)

set(TPLS_meshes mstk metis exodus netcdf)


add_subdirectory(summertest)

# add_executable(MHD MHD.cc PDE_ElectroMagMHD.cc)

# target_link_libraries(MHD ${Amanzi_libs} ${TPLS_meshes} ${Epetra_libs} ${Epetra_interface_libs} ${Teuchos_libs}  
#                           ${Nox_libs} ${Nox_interface_libs} HYPRE superlu_dist superlu parmetis metis
#                           UnitTest++ epetra ${Boost_libs} ${HDF5_libs} ${Amanzi_libs} ${TPLS_meshes} ${Epetra_libs}
#                           ${Epetra_interface_libs} ${Teuchos_libs} ${Nox_libs} ${Nox_interface_libs} HYPRE 
#                           superlu_dist superlu parmetis metis UnitTest++ epetra ${Boost_libs} ${HDF5_libs} ascemio)