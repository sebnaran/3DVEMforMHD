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



#include "PDE_Test.hh"


namespace Amanzi {
namespace Operators {

// Constructor
PDE_Test::PDE_Test(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,Teuchos::ParameterList& plist):
PDE_HelperDiscretization(mesh)
{ 
  //std::cout<<"Entering Constructor"<<std::endl;
  Schema test_schema;
  test_schema.SetBase(AmanziMesh::CELL);
  test_schema.AddItem(AmanziMesh::CELL,WhetStone::DOF_Type::SCALAR,1);
  test_schema.Finalize(mesh); // computes the starting position of the dof ids
  local_op_ = Teuchos::rcp(new Op_Cell_Schema(test_schema,test_schema,mesh) );
  //std::cout<<local_op_->matrices.size()<<std::endl;
  std::cout<<ncells_owned<<std::endl;
  for (int c = 0; c<ncells_owned ; c++) {
    WhetStone::DenseMatrix Acell(1,1);
    Acell(0,0) = 1;
    //std::cout<<Acell<<std::endl;
    //std::cout<<c<<std::endl;
    local_op_->matrices[c] = Acell; 
  }
 // std::cout<<"Created the local matrices"<<std::endl;
  cvs_ = Teuchos::rcp( new CompositeVectorSpace (cvsFromSchema(test_schema, mesh,false)) );
  //std::cout<<"Created the vector space"<<std::endl;
  //std::cout<<global_op_<<std::endl;




  global_op_ = Teuchos::rcp(new Operator_Schema(cvs_, plist, test_schema));
  
  std::cout<<global_op_<<std::endl;
  std::cout<<"Created the global operator"<<std::endl;
  std::cout<<global_op_->rhs()<<std::endl;
  //std::cout<<global_op_->begin()<<::std::endl;
  global_op_->OpPushBack(local_op_);
  //std::cout<<"assigned the local op"<<std::endl;
  global_op_->rhs()->PutScalar(1);
  //std::cout<<"places the number 1 on rhs"<<std::endl;
  std::cout<<"About to leave"<<std::endl;
}


void PDE_Test::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                                    const Teuchos::Ptr<const CompositeVector>& p)
{
  std::cout<<"Matrices updating--not really..."<<std::endl;
}


}  // namespace Operators
}  // namespace Amanzi
