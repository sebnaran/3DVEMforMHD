project(TestSourcesFiles C CXX)

add_executable(ModifiedAmanziTest ModifiedTestMain.cc)

target_link_libraries(ModifiedAmanziTest ${Amanzi_libs} ${TPLS_meshes} ${Epetra_libs} ${Epetra_interface_libs} ${Teuchos_libs}  
                                          ${Nox_libs} ${Nox_interface_libs} HYPRE superlu_dist superlu parmetis metis
                                          UnitTest++ epetra ${Boost_libs} ${HDF5_libs})
add_executable(AssemLinSolv AssembleAndLinearSolve.cc Mesh.cc PDE_Test.cc) #

target_link_libraries(AssemLinSolv ${Amanzi_libs} ${TPLS_meshes} ${Epetra_libs} ${Epetra_interface_libs} ${Teuchos_libs}  
                                           ${Nox_libs} ${Nox_interface_libs} HYPRE superlu_dist superlu parmetis metis
                                           UnitTest++ epetra ${Boost_libs} ${HDF5_libs})

add_executable(TreeOpTest TreeOperatorTest.cc Mesh.cc PDE_Test.cc) #

target_link_libraries(TreeOpTest ${Amanzi_libs} ${TPLS_meshes} ${Epetra_libs} ${Epetra_interface_libs} ${Teuchos_libs}  
                                   ${Nox_libs} ${Nox_interface_libs} HYPRE superlu_dist superlu parmetis metis
                                   UnitTest++ epetra ${Boost_libs} ${HDF5_libs})

add_executable(FirstOrderPoisson FirstOrderPoisson.cc Mesh.cc PDE_FirstOrderPoisson.cc)

target_link_libraries(FirstOrderPoisson ${Amanzi_libs} ${TPLS_meshes} ${Epetra_libs} ${Epetra_interface_libs} ${Teuchos_libs}  
                                 ${Nox_libs} ${Nox_interface_libs} HYPRE superlu_dist superlu parmetis metis
                                 UnitTest++ epetra ${Boost_libs} ${HDF5_libs})
add_executable(SecondOrderPoisson SecondOrderPoisson.cc Mesh.cc PDE_SecondOrderPoisson.cc OutputXDMF.cc)

target_link_libraries(SecondOrderPoisson ${Amanzi_libs} ${TPLS_meshes} ${Epetra_libs} ${Epetra_interface_libs} ${Teuchos_libs}  
                                         ${Nox_libs} ${Nox_interface_libs} HYPRE superlu_dist superlu parmetis metis
                                         UnitTest++ epetra ${Boost_libs} ${HDF5_libs} ${Amanzi_libs} ${TPLS_meshes} ${Epetra_libs}
                                         ${Epetra_interface_libs} ${Teuchos_libs} ${Nox_libs} ${Nox_interface_libs} HYPRE 
                                         superlu_dist superlu parmetis metis UnitTest++ epetra ${Boost_libs} ${HDF5_libs} ascemio)