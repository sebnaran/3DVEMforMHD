#include <UnitTest++.h>

#include "elasticitytest.cc"

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


// Amanzi
#include "MeshFactory.hh"
#include "GMVMesh.hh"
#include "LinearOperatorPCG.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "PDE_Elasticity.hh"

#include "AnalyticElasticity01.hh"
#include "Verification.hh"

//int main(){
//  std::cout<<"Hello"<<std::endl;
//  return 0;
//}


int main(int , const char *[])
{
  return UnitTest::RunAllTests();
}