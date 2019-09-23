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



#include "PDE_ElectroMagMHD.hh"

namespace Amanzi {
namespace Operators {

void PDE_ElectroMagMHD::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                                           const Teuchos::Ptr<const CompositeVector>& p)
{
  std::cout<<"Matrices updating--not really..."<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////////

void PDE_ElectroMagMHD::LocalMagMatrix(){
  Schema magschema;
  magschema.SetBase(AmanziMesh::CELL);
  //The basis used follows the order:
  magschema.AddItem(AmanziMesh::CELL,WhetStone::DOF_Type::SCALAR,3);
  
  //The basis used follows the order:
  magschema.AddItem(AmanziMesh::FACE,WhetStone::DOF_Type::SCALAR,4);
  
  magschema.Finalize(mesh_); // computes the starting position of the dof ids
  //
  mag_local_op_ = Teuchos::rcp(new Op_Cell_Schema(magschema,magschema,mesh_));
  //now we populate the local matrices.
  for (int c = 0; c < ncells_owned; c++){

    unsigned int num_faces = mesh_->cell_get_num_faces(c);
    WhetStone::DenseMatrix Mcell(4*num_faces+3,4*num_faces+3);
    // AmanziMesh::Entity_ID_List faceids;
    // mesh->cell_get_faces(c,&faceids);
    // AmanziGeometry::Point x0,x1,x2,x3;
    // mesh -> node_get_coordinates(nodeids[0],&x0);
    // mesh -> node_get_coordinates(nodeids[1],&x1);
    // mesh -> node_get_coordinates(nodeids[2],&x2);  
    // mesh -> node_get_coordinates(nodeids[3],&x3); 
    //mag_local_op_->matrices[c] = Mcell;
  }
}

}  // namespace Operators
}  // namespace Amanzi
