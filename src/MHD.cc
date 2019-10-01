//Trilinos
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
//Amanzi
#include "VerboseObject_objs.hh"
#include "bilinear_form_registration.hh"
#include "MeshFactory.hh"
#include "LinearOperatorPCG.hh"
#include "TreeOperator.hh"
#include "OutputXDMF.hh"
#include "Polynomial.hh"
//Proper
#include "PDE_ElectroMagMHD.hh"
#include "PolynomialByCoordinates.hh"
#include "InitialMagneticField.hh"
int main(int argc, char *argv[]){
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;
  
  Teuchos::GlobalMPISession mpiSession(&argc, &argv); 
  auto comm = Amanzi::getDefaultComm();
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
 // Generate a hex mesh (3D)
  //
  // Generates a structured mesh covering [x0,x1] X [y0,y1] X [z0,z1] with
  // (nx, ny, nz) cells.
  // Teuchos::RCP<Mesh> create(const double x0, const double y0, const double z0,
  //                           const double x1, const double y1, const double z1,
  //                           const int nx, const int ny, const int nz,
  //                           const bool request_faces=true,
  //                           const bool request_edges=false);
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(-1,-1,-1,1,1,1,3,3,3);
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nnodes = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);
  
  Teuchos::RCP<PDE_ElectroMagMHD> electro_op = Teuchos::rcp(new PDE_ElectroMagMHD(mesh));
  electro_op->LocalMagMatrix();
  InitialMagneticField Bz;
  // WhetStone::NumericalIntegration numi(mesh);

  // PolynomialByCoordinates pol(3,2);
  // std::cout<<"The degree is"<<pol.GetDeg()<<"and the dim is"<<pol.GetDim()<<std::endl;
  // auto poly = pol.GetPol1();

  // WhetStone::Polynomial pol1(3,2);
  // pol1(1,1) = 5;
  // pol1(2,3) = 6;
  // pol1(1,3) = 6;
  // std::cout<<pol1<<std::endl;
  // PolynomialByCoordinates Polys1(3,2,pol1,pol1,pol1);
  // PolynomialByCoordinates Polys2(3,2,pol1,pol1,pol1);
  // Amanzi::WhetStone::Polynomial polyf = Polys1.dot(Polys2);
  // std::cout<<polyf<<std::endl;
  // WhetStone::Polynomial poly(3,1);
  // double val;
  // poly(1, 0) = 2.0;
  // poly(1, 1) = 3.0;
  // poly(1, 2) = 4.0;
  // int cell(0);
  // std::cout<<poly<<std::endl;
  // std::vector<const WhetStone::WhetStoneFunction*> polys(1);
  // polys[0] = &poly;
  // val = numi.IntegrateFunctionsTriangulatedCell(cell, polys, poly.order());
  // std::cout<<"The integral over the cell is "<<val<<std::endl;
  // poly.set_origin(mesh->cell_centroid(cell));
  // std::cout<<poly<<std::endl;
  // val = numi.IntegrateFunctionsTriangulatedCell(cell, polys, poly.order());
  // std::cout<<"The integral over the cell is "<<val<<std::endl;
   //Next we will define the boundary conditions
  // Teuchos::RCP<BCs> bcv = Teuchos::rcp(new BCs(mesh, AmanziMesh::NODE, WhetStone::DOF_Type::SCALAR));
  // std::vector<int>& bcv_model = bcv->bc_model();
  // std::vector<double>& bcv_value = bcv->bc_value();
  
  // Point xv(2);
  // for (int v = 0; v < nnodes_wghost; ++v) {
  //  mesh->node_get_coordinates(v, &xv);
  // //We will identify which points lie in the boundary
  // if (fabs(xv[0]+1) < 1e-6 || fabs(xv[0] - 1.0) < 1e-6 ||
  //     fabs(xv[1]+1) < 1e-6 || fabs(xv[1] - 1.0) < 1e-6) {
  //     bcv_model[v] = OPERATOR_BC_DIRICHLET;
  //     bcv_value[v] = 25;
  //   }}

  //Teuchos::RCP<PDE_SecondOrderPoisson> electro_op = Teuchos::rcp(new PDE_SecondOrderPoisson(mesh));
  //The global stiffness operator
  // Teuchos::RCP<Operator> global_op = op_poisson->global_operator();
  // const CompositeVectorSpace & cvs = *op_poisson->GetCVS();
  // //Before applying the boundary conditions we must populate the RHS
  // CompositeVector source(cvs);
  // Epetra_MultiVector& src = *source.ViewComponent("node");
  // for (int v = 0; v < nnodes; v++) {
  //   mesh->node_get_coordinates(v, &xv);
  //   src[0][v] = 1;
  // }
  // global_op->UpdateRHS(source, true);
  // op_poisson->SetBCs(bcv, bcv);   
  // //Global assembly
  // op_poisson->ApplyBCs(true, true, true);
  // global_op->SymbolicAssembleMatrix();
  // global_op->AssembleMatrix();
  // //auto S = *global_op->A();
  // //std::cout<<S<<std::endl;
  // //Initialize the solution
  // CompositeVector solution(cvs);
  // solution.PutScalar(0.0);

  // std::string xmlFileName = "test/operator_elasticity.xml";
  // Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  // Teuchos::ParameterList plist = xmlreader.getParameters();
  // Teuchos::ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
  // global_op->InitializePreconditioner(slist);
  // global_op->UpdatePreconditioner();

  // Teuchos::ParameterList lop_list = plist.sublist("solvers")
  //                                        .sublist("PCG").sublist("pcg parameters");
  // AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
  //     pcg(global_op, global_op);
  // pcg.Init(lop_list);

  // CompositeVector& rhs = *global_op->rhs();
  // int ierr = pcg.ApplyInverse(rhs, solution);
  // std::cout<<"sol is ="<<std::endl;
  // solution.Print(std::cout);
  // std::cout<<"rhs is"<<std::endl;
  // rhs.Print(std::cout);
  // std::cout<<"The solution is"<<std::endl;
  // solution.Print(std::cout);
  // Teuchos::ParameterList iolist;
  // iolist.get<std::string>("file name base", "plot");
  // OutputXDMF io(iolist, mesh, true, false);
  // io.InitializeCycle(0,0);
  // const Epetra_MultiVector & sol = *solution.ViewComponent("node",false);
  // std::cout<<"The epetra version of sol is"<<std::endl;
  // std::cout<<sol<<std::endl;
  // io.WriteVector(*sol(0), "solution", AmanziMesh::NODE);
  // io.FinalizeCycle();
}