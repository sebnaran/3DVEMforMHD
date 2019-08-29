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



#include "PDE_FirstOrderPoisson.hh"


namespace Amanzi {
namespace Operators {

// Constructor
PDE_FirstOrderPoisson::PDE_FirstOrderPoisson(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)://,Teuchos::ParameterList& plist):
PDE_HelperDiscretization(mesh)
{
  //std::cout<<"Entering Constructor"<<std::endl;
  Schema u_schema;
  Schema p_schema;
  //The assembly should run over the cells thus its base are the cells
  p_schema.SetBase(AmanziMesh::CELL);
  //The pressure DOFs are cell-based and scalars.
  p_schema.AddItem(AmanziMesh::CELL,WhetStone::DOF_Type::SCALAR,1);
  p_schema.Finalize(mesh); // computes the starting position of the dof ids
  //
  //The assembly should run over the cells thus its base are the cells
  u_schema.SetBase(AmanziMesh::CELL);
  //The velocity DOFs are face-based and scalars
  u_schema.AddItem(AmanziMesh::FACE,WhetStone::DOF_Type::SCALAR,1);
  u_schema.Finalize(mesh);
  //
  //Next we will create a container for the local matrices for the rectangular
  //divergence and the square mass matrices.
   mass_local_op_ = Teuchos::rcp(new Op_Cell_Schema(u_schema,u_schema,mesh));
  //Now we populate the local matrices
   for (int c = 0; c < ncells_owned ; c++) {
     WhetStone::DenseMatrix Mcell(2,2);
     Mcell(0,0) = 1.00, Mcell(0,1) = 1.00;
     Mcell(1,0) = 1.00, Mcell(1,1) = 1.00;
     mass_local_op_->matrices[c] = Mcell;
    }
  // Op_Cell_Schema(row,column,mesh)we will put the velocities in the columns
  // pressures on the rows so the matrix we construct will aproximate:
  // int_P q div u = q^T S u
  //The constructor for Op_Cell_Schema is given below
  // Op_Cell_Schema(const Schema& schema_row, const Schema& schema_col,
  //                const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
  //     Op(schema_row, schema_col, mesh)
    div_local_op_ = Teuchos::rcp(new Op_Cell_Schema(p_schema,u_schema,mesh));
  //Now we populate the local matrices
     for (int c = 0; c < ncells_owned ; c++){
     WhetStone::DenseMatrix Scell(1,4);
     Scell(0,0) = 1, Scell(0,1) = 1, Scell(0,2) = 1, Scell(0,3) = 1;
     div_local_op_->matrices[c] = Scell;
    }
    //Now we can define the velocity and pressure spaces
    p_cvs_ = Teuchos::rcp( new CompositeVectorSpace (cvsFromSchema(p_schema, mesh)) );
    u_cvs_ = Teuchos::rcp( new CompositeVectorSpace (cvsFromSchema(u_schema, mesh)) );
    //The constructor for a global operator requires a parameter list so we will create a dummy one
    Teuchos::ParameterList plist = Teuchos::ParameterList();
    //This line creates a global operator for the mass matrix
    mass_global_op_ = Teuchos::rcp(new Operator_Schema(u_cvs_, plist, u_schema));
    //This line assigns the corresponding container of local matrices.
    mass_global_op_->OpPushBack(mass_local_op_);
    //The constructor for operator scheme looks like this:
    // general rectangular operator
    // Operator_Schema(const Teuchos::RCP<const CompositeVectorSpace>& cvs_row,
    //               const Teuchos::RCP<const CompositeVectorSpace>& cvs_col,
    //               Teuchos::ParameterList& plist,
    //               const Schema& schema_row,
    //               const Schema& schema_col)

    div_global_op_ = Teuchos::rcp(new Operator_Schema(p_cvs_,u_cvs_,plist,p_schema,u_schema));
    div_global_op_->OpPushBack(div_local_op_);
}


void PDE_FirstOrderPoisson::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                                           const Teuchos::Ptr<const CompositeVector>& p)
{
  std::cout<<"Matrices updating--not really..."<<std::endl;
}


}  // namespace Operators
}  // namespace Amanzi
