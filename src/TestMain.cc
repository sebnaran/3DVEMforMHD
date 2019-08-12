#include <UnitTest++.h>

#include "elasticitytest.cc"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "mpi.h"

// TPLs
#include "EpetraExt_RowMatrixOut.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_GlobalMPISession.hpp"

// Amanzi
#include "MeshFactory.hh"
#include "GMVMesh.hh"
#include "LinearOperatorPCG.hh"
#include "Tensor.hh"
#include "VerboseObject_objs.hh"

// Amanzi::Operators
#include "PDE_Elasticity.hh"

#include "AnalyticElasticity01.hh"
#include "Verification.hh"

//int main(){
//  std::cout<<"Hello"<<std::endl;
//  return 0;
//}

TEST(Sanity) 
{
   CHECK_EQUAL(1, 1);
}


int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return UnitTest::RunAllTests();
}
