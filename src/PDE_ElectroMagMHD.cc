#include <vector>

// TPLs
#include "Epetra_Vector.h"

// Amanzi
//#include "BilinearFormFactory.hh"
//#include "bilinear_form_registration.hh"
#include "errors.hh"
#include "MatrixFE.hh"
#include "PreconditionerFactory.hh"
#include "WhetStoneDefs.hh"
#include "MeshDefs.hh"
#include "Polynomial.hh"
#include "NumericalIntegration.hh"
// Amanzi::Operators
#include "Op.hh"
#include "Op_Cell_Schema.hh"
#include "OperatorDefs.hh"
#include "Operator_Schema.hh"
#include "PDE_Elasticity.hh"
#include "Schema.hh"



#include "PDE_ElectroMagMHD.hh"
#include "PolynomialByCoordinates.hh"
namespace Amanzi {
namespace Operators {

PDE_ElectroMagMHD::PDE_ElectroMagMHD(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh):
  PDE_HelperDiscretization(mesh){
  //Creating the magnetic schema
  magschema_.SetBase(AmanziMesh::CELL);
  //The cell degrees of freedom are moments agaist these
  //pg_0 = (0,0,yP) pg_1 = (zP,0,0) pg_2=(0,xP,0)
  magschema_.AddItem(AmanziMesh::CELL,WhetStone::DOF_Type::SCALAR,3);
  
  //The basis used follows the order:
  magschema_.AddItem(AmanziMesh::FACE,WhetStone::DOF_Type::SCALAR,4);
  magschema_.Finalize(mesh_); // computes the starting position of the dof ids
  magcvs_ = Teuchos::rcp( new CompositeVectorSpace (cvsFromSchema(magschema_, mesh_,false)) );
  mag_local_op_ = Teuchos::rcp(new Op_Cell_Schema(magschema_,magschema_,mesh_));
  }


void PDE_ElectroMagMHD::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                                           const Teuchos::Ptr<const CompositeVector>& p)
{
  std::cout<<"Matrices updating--not really..."<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////////

void PDE_ElectroMagMHD::LocalMagMatrix(){
  //now we populate the local matrices.
  //for (int c = 0; c < ncells_owned; c++){
  for (int c = 0; c < 1; c++){
    WhetStone::Polynomial zero(3,0);
    zero.set_origin(mesh_->cell_centroid(c));
    zero(0,0) = 0;

    WhetStone::Polynomial one(3,0);
    one.set_origin(mesh_->cell_centroid(c));
    one(0,0) = 1;

    WhetStone::Polynomial x(3,1);
    x.set_origin(mesh_->cell_centroid(c));
    x(1,0) = 1;

    WhetStone::Polynomial y(3,1);
    y.set_origin(mesh_->cell_centroid(c));
    y(1,1) = 1;

    WhetStone::Polynomial z(3,1);
    z.set_origin(mesh_->cell_centroid(c));
    z(1,2) = 1;
    
        //The polynomial basis we use is
    //q_0 = (0,0,yP), q_1 = (zP,0,0) q_2 = (0,xP,0)
    //q_3 = (1,0,0),  q_4 = (xP,0,0) q_5 = (0,1,0)
    //q_6 = (0,yP,0), q_7 = (0,0,1)  q_8 = (0,0,zP)
    //q_9 = (0,zP,yP), q_10 = (zP,0,xP), q_11 = (yP,xP,0)
    PolynomialByCoordinates q0(3,1,zero,zero,y);
    PolynomialByCoordinates q1(3,1,z,zero,zero);
    PolynomialByCoordinates q2(3,1,zero,x,zero);
    PolynomialByCoordinates q3(3,0,one,zero,zero);
    PolynomialByCoordinates q4(3,1,x,zero,zero);
    PolynomialByCoordinates q5(3,0,zero,one,zero);
    PolynomialByCoordinates q6(3,1,zero,y,zero);
    PolynomialByCoordinates q7(3,0,zero,zero,one);
    PolynomialByCoordinates q8(3,1,zero,zero,z);
    PolynomialByCoordinates q9(3,1,zero,z,y);
    PolynomialByCoordinates q10(3,1,z,zero,x);
    PolynomialByCoordinates q11(3,1,y,x,zero);

    std::vector<PolynomialByCoordinates*> basis{&q0,&q1,&q2,&q3,&q4,&q5,&q6,&q7,&q8,&q9,&q10,&q11};
    WhetStone::NumericalIntegration numi(mesh_); 
    std::vector<const WhetStone::WhetStoneFunction*> polys(1);
    unsigned int num_faces = mesh_->cell_get_num_faces(c);
    unsigned int num_local_dof = 4*num_faces+3;
    WhetStone::DenseMatrix D(num_local_dof,12);
    WhetStone::DenseMatrix G(12,12);
    WhetStone::Polynomial poly;
    for(int i = 0;i < 12;i++){
      for(int j = 0;j <12;j++){
        //std::cout<<"("<<i<<","<<j<<")"<<std::endl;
        //std::cout<<basis[i]->dot(*basis[j])<<std::endl;
        poly = basis[i]->dot(*basis[j]);
        polys[0] = &poly;
        G(i,j) = numi.IntegrateFunctionsTriangulatedCell(c, polys, poly.order());
        if (i<3){
          D(i,j) = G(i,j);  
        }
      }
    }
    AmanziMesh::Entity_ID_List faceids;
    mesh_->cell_get_faces(c,&faceids);

  //  for(int i = 3;i < num_local_dof;i++){
  //     for(int j = 0;j <12;j++){
  //         D(i,j) = 1;
  //         std::cout<<faceids[i]<<std::endl;
  //     }
  //  }




    // //row 0
    // poly(2,3) = 1;
    // poly.set_origin(mesh_->cell_centroid(c));
    // std::vector<const WhetStone::WhetStoneFunction*> polys(1);
    // polys[0] = &poly;
    // G(0,0) = numi.IntegrateFunctionsTriangulatedCell(c, polys, poly.order());
    // // (0,1)-(0,7)
    // G(0,1) = 0, G(0,2) = 0, G(0,3) = 0;
    // G(0,4) = 0, G(0,5) = 0, G(0,6) = 0;
    // G(0,7) = 0;
    // G(1,0) = 0, G(2,0) = 0, G(3,0) = 0;
    // G(4,0) = 0, G(5,0) = 0, G(6,0) = 0;
    // G(7,0) = 0;
    // //
    // poly.Reshape(3,2,true);
    // poly(2,4) = 1; 
 
    // G(0,8) = numi.IntegrateFunctionsTriangulatedCell(c, polys, poly.order());
    // G(8,0) = G(0,8);
    
    // G(0,9) = G(0,0);
    // G(9,0) = G(0,9);

    // poly.Reshape(3,2,true);
    // poly(2,1) = 1;

    // G(0,10) = numi.IntegrateFunctionsTriangulatedCell(c, polys, poly.order()); 
    // G(10,0) = G(0,10);
    // G(0,11) = 0, G(11,0) = 0;
    
    // //Row 1
    //   //The polynomial basis we use is
    // //q_0 = (0,0,yP), q_1 = (zP,0,0) q_2 = (0,xP,0)
    // //q_3 = (1,0,0),  q_4 = (xP,0,0) q_5 = (0,1,0)
    // //q_6 = (0,yP,0), q_7 = (0,0,1)  q_8 = (0,0,zP)
    // //q_9 = (0,zP,yP), q_10 = (zP,0,xP), q_11 = (yP,xP,0)
    // poly.Reshape(3,2,true);
    // poly(2,5) = 1;
    // G(1,1) = numi.IntegrateFunctionsTriangulatedCell(c, polys, poly.order());
    
    // G(1,2) = 0, G(2,1) = 0;
    // G(1,3) = 0, G(1,3) = 0;

    // poly.Reshape(3,2,true);
    // poly(2,2) = 1;
    // G(1,4) = numi.IntegrateFunctionsTriangulatedCell(c, polys, poly.order());
    // G(4,1) = G(1,4);

    // G(1,5) = 0, G(1,6) = 0, G(1,7) = 0, G(1,8) = 0, G(1,9) = 0;
    // G(5,1) = 0, G(6,1) = 0, G(7,1) = 0, G(8,1) = 0, G(9,1) = 0;
    
    // G(1,10) = G(1,1), G(10,1) = G(1,10);

    // poly.Reshape(3,2,true);
    // poly(2,4) = 1;
    // G(1,11) = numi.IntegrateFunctionsTriangulatedCell(c, polys, poly.order());
    // G(11,1) = G(1,11);
    
    // //row 2
    // //The polynomial basis we use is
    // //q_0 = (0,0,yP), q_1 = (zP,0,0) q_2 = (0,xP,0)
    // //q_3 = (1,0,0),  q_4 = (xP,0,0) q_5 = (0,1,0)
    // //q_6 = (0,yP,0), q_7 = (0,0,1)  q_8 = (0,0,zP)
    // //q_9 = (0,zP,yP), q_10 = (zP,0,xP), q_11 = (yP,xP,0)
    
    // poly.Reshape(3,2,true);
    // poly(2,0) = 1;
    // std::cout<<poly<<std::endl; 
    // G(2,2) = numi.IntegrateFunctionsTriangulatedCell(c, polys, poly.order());
    // G(11,1) = G(1,11);

    // //std::cout<<G(0,8)<<std::endl;
    // // AmanziMesh::Entity_ID_List faceids;
    // // mesh->cell_get_faces(c,&faceids);
    // // AmanziGeometry::Point x0,x1,x2,x3;
    // // mesh -> node_get_coordinates(nodeids[0],&x0);
    // // mesh -> node_get_coordinates(nodeids[1],&x1);
    // // mesh -> node_get_coordinates(nodeids[2],&x2);  
    // // mesh -> node_get_coordinates(nodeids[3],&x3); 
    // //mag_local_op_->matrices[c] = Mcell;
    // unsigned int num_faces = mesh_->cell_get_num_faces(c);
    // WhetStone::DenseMatrix Mcell(4*num_faces+3,4*num_faces+3);
  }
 
  //void PDE_ElectroMagMHD::ContructLocalMatrix(){
    
  //}
   
  
}

CompositeVector PDE_ElectroMagMHD::MagneticDOFs(WhetStone::WhetStoneFunction* Bx,
                                                WhetStone::WhetStoneFunction* By,
                                                WhetStone::WhetStoneFunction* Bz){
CompositeVector Bh(*magcvs_);
std::cout<<Bh.NumComponents()<<std::endl;
return Bh;
}



}  // namespace Operators
}  // namespace Amanzi
