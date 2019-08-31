#include <vector>

// TPLs
#include "Epetra_Vector.h"

// Amanzi
#include "BilinearFormFactory.hh"
#include "errors.hh"
#include "MatrixFE.hh"
#include "PreconditionerFactory.hh"
#include "WhetStoneDefs.hh"
#include "MeshDefs.hh"
// Amanzi::Operators
#include "Op.hh"
#include "Op_Cell_Schema.hh"
#include "OperatorDefs.hh"
#include "Operator_Schema.hh"
#include "PDE_Elasticity.hh"
#include "Schema.hh"



#include "PDE_SecondOrderPoisson.hh"


namespace Amanzi {
namespace Operators {

// Constructor
PDE_SecondOrderPoisson::PDE_SecondOrderPoisson(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)://,Teuchos::ParameterList& plist):
PDE_HelperDiscretization(mesh)
{
  //std::cout<<"Entering Constructor"<<std::endl;
  Schema p_schema;
  //The assembly should run over the cells thus its base are the cells
  p_schema.SetBase(AmanziMesh::CELL);
  //The pressure DOFs are cell-based and scalars.
  p_schema.AddItem(AmanziMesh::NODE,WhetStone::DOF_Type::SCALAR,1);
  p_schema.Finalize(mesh); // computes the starting position of the dof ids
  //
  //Next we will create a container for the local matrices for the rectangular
  //divergence and the square mass matrices.
   local_op_ = Teuchos::rcp(new Op_Cell_Schema(p_schema,p_schema,mesh));
  //Now we populate the local matrices
   for (int c = 0; c < ncells_owned ; c++) {
     double area = mesh->cell_area(c);  
     WhetStone::DenseMatrix Mcell(4,4);
     
     Mcell(0,0) = 1.00, Mcell(0,1) = 1.00, Mcell(0,2) = 1, Mcell(0,3) = 1;
     Mcell(1,0) = 1.00, Mcell(1,1) = 1.00, Mcell(1,2) = 1, Mcell(1,3) = 1;
     Mcell(2,0) = 1.00, Mcell(2,1) = 1.00, Mcell(2,2) = 1, Mcell(2,3) = 1;
     Mcell(3,0) = 1.00, Mcell(3,1) = 1.00, Mcell(3,2) = 1, Mcell(3,3) = 1;
     local_op_->matrices[c] = Mcell;
    }
    //Now we can define the velocity and pressure spaces
    cvs_ = Teuchos::rcp( new CompositeVectorSpace (cvsFromSchema(p_schema, mesh)) );
    //The constructor for a global operator requires a parameter list so we will create a dummy one
    Teuchos::ParameterList plist = Teuchos::ParameterList();
    //This line creates a global operator for the mass matrix
    global_op_ = Teuchos::rcp(new Operator_Schema(cvs_, plist, p_schema));
    //This line assigns the corresponding container of local matrices.
    global_op_->OpPushBack(local_op_);
}


void PDE_SecondOrderPoisson::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                                           const Teuchos::Ptr<const CompositeVector>& p)
{
  std::cout<<"Matrices updating--not really..."<<std::endl;
}


}  // namespace Operators
}  // namespace Amanzi
