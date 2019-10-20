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
  p_schema.set_base(AmanziMesh::CELL);
  //The pressure DOFs are cell-based and scalars.
  p_schema.AddItem(AmanziMesh::NODE,WhetStone::DOF_Type::SCALAR,1);
  p_schema.Finalize(mesh); // computes the starting position of the dof ids
  //
  //Next we will create a container for the local matrices 
   local_op_ = Teuchos::rcp(new Op_Cell_Schema(p_schema,p_schema,mesh));
  //These are the structures that we require to build the local matrices
   AmanziMesh::Entity_ID_List edgeids;
   AmanziMesh::Entity_ID Node0,Node1,Node2,Node3;
   AmanziMesh::Entity_ID_List nodeids;

   AmanziGeometry::Point v0,v1,v2,v3;
   AmanziGeometry::Point n0,n1,n2,n3;
   AmanziGeometry::Point vP;

   WhetStone::DenseMatrix N(4,3);
   WhetStone::DenseMatrix R(4,3);
   WhetStone::DenseMatrix TR(3,4);
   WhetStone::DenseMatrix Mcell(4,4);
   WhetStone::DenseMatrix Scell(4,4);
   WhetStone::DenseMatrix temp1(3,3);
   WhetStone::DenseMatrix temp2(4,3);

   double lambda;
   double l0,l1,l2,l3;
   //Now we will populate these matrices
   for (int c = 0; c < ncells_owned ; c++) {
    
    vP = mesh_->cell_centroid(c);
    mesh -> cell_get_edges(c, &edgeids);
    
    mesh_ -> edge_get_nodes(edgeids[0],&Node0,&Node1);
    mesh_ -> edge_get_nodes(edgeids[1],&Node1,&Node2);
    mesh_ -> edge_get_nodes(edgeids[2],&Node2,&Node3);
   
    mesh_ -> node_get_coordinates(Node0,&v0);
    mesh_ -> node_get_coordinates(Node1,&v1);
    mesh_ -> node_get_coordinates(Node2,&v2);  
    mesh_ -> node_get_coordinates(Node3,&v3); 

     N(0,0) = 1, N(0,1) = v0[0]-vP[0], N(0,2) = v0[1]-vP[1];  
     N(1,0) = 1, N(1,1) = v1[0]-vP[0], N(1,2) = v1[1]-vP[1];
     N(2,0) = 1, N(2,1) = v2[0]-vP[0], N(2,2) = v2[1]-vP[1];
     N(3,0) = 1, N(3,1) = v3[0]-vP[0], N(3,2) = v3[1]-vP[1];
     
     l0 = mesh_ -> edge_length(edgeids[0]);
     l1 = mesh_ -> edge_length(edgeids[1]);
     l2 = mesh_ -> edge_length(edgeids[2]);
     l3 = mesh_ -> edge_length(edgeids[3]);

     n0 = mesh_ -> face_normal(edgeids[0],false,c);
     n1 = mesh_ -> face_normal(edgeids[1],false,c);
     n2 = mesh_ -> face_normal(edgeids[2],false,c);
     n3 = mesh_ -> face_normal(edgeids[3],false,c);

     R(0,0) = 0, R(0,1) = l3*n3[0]+l0*n0[0], R(0,2) = l3*n3[1]+l0*n0[1];
     R(1,0) = 0, R(1,1) = l0*n0[0]+l1*n1[0], R(1,2) = l0*n0[1]+l1*n1[1];
     R(2,0) = 0, R(2,1) = l1*n1[0]+l2*n2[0], R(2,2) = l1*n1[1]+l2*n2[1];
     R(3,0) = 0, R(3,1) = l2*n2[0]+l3*n3[0], R(3,2) = l2*n2[1]+l3*n3[1];

    //Mcell = R(R^TN)^t R^T+[tr(R(R^TN)^t R^T)]*(Id-N(N^TN)^{-1}N^T)/2

    int k1 = temp1.Multiply(R,N,true);
    int k2 = temp1.InverseMoorePenrose();
    int k3 = temp2.Multiply(R,temp1,false);
    //std::cout<<"R="<<R<<std::endl;
    TR.Transpose(R);
    //std::cout<<"TR="<<TR<<std::endl;
    int k4 = Mcell.Multiply(temp2,TR,false);
    lambda = temp1.Trace();
    
    int d1 = temp1.Multiply(N,N,true);
    int d2 = temp1.InverseSPD();
    int d3 = temp2.Multiply(N,temp1,false);
    //std::cout<<"N="<<N<<std::endl;
    TR.Transpose(N);
    //std::cout<<"TR="<<TR<<std::endl;
    int d4 = Scell.Multiply(temp2,TR,false);
    for(int i=0;i<4;i++){
      for(int j=0;j<4;j++){
         if(i==j){
           Scell(i,j) = lambda*(1-Scell(i,j))/2;
         }
         else{
           Scell(i,j) = -lambda*(Scell(i,j))/2;
         }
      }
    }
    local_op_->matrices[c] = Mcell+Scell;
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
