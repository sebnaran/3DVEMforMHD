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
  //Next we will create a container for the local matrices 
   local_op_ = Teuchos::rcp(new Op_Cell_Schema(p_schema,p_schema,mesh));
  //Now we populate the local matrices
   AmanziMesh::Entity_ID_List edgeids;
   for (int c = 0; c < ncells_owned ; c++) {
    // AmanziGeometry::Point vP = mesh_->cell_centroid(c);
    //std::cout<<mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED)<<std::endl;
    mesh -> cell_get_edges(c, &edgeids);
    //  AmanziMesh::Entity_ID Node0,Node1,Node2,Node3;
    //  mesh_ -> edge_get_nodes(edgeids[0],&Node0,&Node1);
    //  mesh_ -> edge_get_nodes(edgeids[1],&Node1,&Node2);
    //  mesh_ -> edge_get_nodes(edgeids[2],&Node2,&Node3);
     //AmanziMesh::Entity_ID_List nodeids;
    // mesh->cell_get_nodes(c,&nodeids);
    //  AmanziGeometry::Point v0,v1,v2,v3;
    //  mesh_ -> node_get_coordinates(Node0,&v0);
    //  mesh_ -> node_get_coordinates(Node1,&v1);
    //  mesh_ -> node_get_coordinates(Node2,&v2);  
    //  mesh_ -> node_get_coordinates(Node3,&v3); 

    //  WhetStone::DenseMatrix N(4,3);
    //  N(0,0) = 1, N(0,1) = v0[0]-vP[0], N(0,2) = v0[1]-vP[1];  
    //  N(1,0) = 1, N(1,1) = v1[0]-vP[0], N(1,2) = v1[1]-vP[1];
    //  N(2,0) = 1, N(2,1) = v2[0]-vP[0], N(2,2) = v2[1]-vP[1];
    //  N(3,0) = 1, N(3,1) = v3[0]-vP[0], N(3,2) = v3[1]-vP[1];
     
    //  double l0 = mesh_ -> edge_length(edgeids[0]);
    //  double l1 = mesh_ -> edge_length(edgeids[1]);
    //  double l2 = mesh_ -> edge_length(edgeids[2]);
    //  double l3 = mesh_ -> edge_length(edgeids[3]);

    //  AmanziGeometry::Point n0,n1,n2,n3;
    //  n0 = mesh_ -> face_normal(edgeids[0],false,c);
    //  n1 = mesh_ -> face_normal(edgeids[1],false,c);
    //  n2 = mesh_ -> face_normal(edgeids[2],false,c);
    //  n3 = mesh_ -> face_normal(edgeids[3],false,c);
    //  std::cout<<n0<<std::endl;
    //  WhetStone::DenseMatrix R(4,3);

     //R(0,0) = 0, R(0,1) = , R(0,2)= 
     WhetStone::DenseMatrix Mcell(4,4);
     Mcell(0,0) = 1.00, Mcell(0,1) = 1.00, Mcell(0,2) = 1, Mcell(0,3) = 1;
     Mcell(1,0) = 1.00, Mcell(1,1) = 1.00, Mcell(1,2) = 1, Mcell(1,3) = 1;
     Mcell(2,0) = 1.00, Mcell(2,1) = 1.00, Mcell(2,2) = 1, Mcell(2,3) = 1;
     Mcell(3,0) = 1.00, Mcell(3,1) = 1.00, Mcell(3,2) = 1, Mcell(3,3) = 1;
     //Mcell(0,0) = 2*(pow(area,2.0)/16), Mcell(0,1) = 0;
     local_op_->matrices[c] = Mcell;
    }
    //Now we can define the velocity and pressure spaces
    cvs_ = Teuchos::rcp( new CompositeVectorSpace (cvsFromSchema(p_schema, mesh,false)) );
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
