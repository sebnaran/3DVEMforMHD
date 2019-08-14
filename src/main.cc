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
    std::cout<<"Before argc = "<<argc<<std::endl;
    std::cout<<"Before argv = "<<argv<<std::endl;
    std::cout<<"Before *argv = "<<*argv[0]<<std::endl;
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);
    std::cout<<"After GlobalMpi session argc = "<<argc<<std::endl;
    std::cout<<"After GlobalMpi session argv = "<<argv<<std::endl;
    std::cout<<"After GlobalMpi session *argv = "<<*argv[0]<<std::endl;
    auto comm = Amanzi::getDefaultComm();
    std::cout<<"After Comm argc = "<<argc<<std::endl;
    std::cout<<"After Comm argv = "<<argv<<std::endl;
    std::cout<<"After Comm *argv = "<<*argv[0]<<std::endl;
    std::cout<<"Comm="<<comm<<std::endl;
    int MyPID = comm->MyPID();
    std::cout<<"MyPID="<<MyPID<<std::endl;

    std::string xmlFileName = "test/operator_elasticity.xml";
    Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
    Teuchos::ParameterList plist = xmlreader.getParameters();
    Teuchos::ParameterList op_list = plist.sublist("PK operator")
                                        .sublist("elasticity operator");
    
    std::cout<< "plist="<<plist<<std::endl;   
    std::cout<< "op_list="<<op_list<<std::endl; 
    
    return 0;
}