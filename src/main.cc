#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "EpetraExt_RowMatrixOut.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "GMVMesh.hh"
#include "LinearOperatorPCG.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "PDE_Elasticity.hh"
#include "Verification.hh"
#include "AnalyticElasticity01.hh"

int main(int argc, char *argv[]){
//     Teuchos::GlobalMPISession mpiSession(&argc, &argv);
     auto comm = Amanzi::getDefaultComm();
//     int MyPID = comm->MyPID();
//   if (MyPID == 0) std::cout << "\nTest: 2D elasticity: exactness test" << std::endl;

    // std::string xmlFileName = "test/operator_elasticity.xml";
    // Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
    // Teuchos::ParameterList plist = xmlreader.getParameters();
    // Teuchos::ParameterList interlist = plist.sublist("PK operator");
    // Teuchos::ParameterList op_list = plist.sublist("PK operator")
    //                                     .sublist("elasticity operator");
    // std::cout<<"plist is "<<std::endl<<plist<<std::endl;
    // std::cout<<"interlist is "<<std::endl<<interlist<<std::endl;
    // std::cout<<"op_list is "<<std::endl<<op_list<<std::endl;
  using namespace Amanzi::AmanziMesh;
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 4, 5);

  // -- general information about mesh
//   int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
//   int nnodes = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
//   int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
//   int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
//   int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);


    return 0;
}