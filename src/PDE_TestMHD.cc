/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

// TPLs
#include "Epetra_Vector.h"

// Amanzi
#include "BilinearFormFactory.hh"
#include "errors.hh"
#include "MatrixFE.hh"
#include "PreconditionerFactory.hh"
#include "WhetStoneDefs.hh"

// Amanzi::Operators
#include "Op.hh"
#include "Op_Cell_Schema.hh"
#include "OperatorDefs.hh"
#include "Operator_Schema.hh"
#include "PDE_Elasticity.hh"
#include "PDE_TestMHD.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

void Test_MHD::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                                    const Teuchos::Ptr<const CompositeVector>& p)
{
  std::cout<<"Matrices updating--not really..."<<std::endl;
}


}  // namespace Operators
}  // namespace Amanzi
