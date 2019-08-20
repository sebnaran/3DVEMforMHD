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
#include "Schema.hh"

#include "PDE_Test.hh"


namespace Amanzi {
namespace Operators {

// Constructor
PDE_Test::PDE_Test(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh):
PDE_HelperDiscretization(mesh){}


void PDE_Test::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                                    const Teuchos::Ptr<const CompositeVector>& p)
{
  std::cout<<"Matrices updating--not really..."<<std::endl;
}


}  // namespace Operators
}  // namespace Amanzi
