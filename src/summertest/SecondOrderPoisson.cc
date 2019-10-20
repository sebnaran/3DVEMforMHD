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

//Proper
#include "PDE_SecondOrderPoisson.hh"

int main(int argc, char *argv[]){
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;
  
  Teuchos::GlobalMPISession mpiSession(&argc, &argv); 
  auto comm = Amanzi::getDefaultComm();
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
  // Generates a structured mesh covering [x0,x1] X [y0,y1] with (nx, ny)
  // cells.
  //  Teuchos::RCP<Mesh> create(const double x0, const double y0,
  //                           const double x1, const double y1,
  //                           const int nx, const int ny,
  //                           const bool request_faces=true,
  //                           const bool request_edges=false);
  int n = 20;
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(-1,-1,1,1,n,n,true,true);
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nnodes = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);
  
   //Next we will define the boundary conditions
  Teuchos::RCP<BCs> bcv = Teuchos::rcp(new BCs(mesh, AmanziMesh::NODE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bcv_model = bcv->bc_model();
  std::vector<double>& bcv_value = bcv->bc_value();
  
  Point xv(2);
  for (int v = 0; v < nnodes_wghost; ++v) {
   mesh->node_get_coordinates(v, &xv);
  //We will identify which points lie in the boundary
  if (fabs(xv[0]+1) < 1e-6 || fabs(xv[0] - 1.0) < 1e-6 ||
      fabs(xv[1]+1) < 1e-6 || fabs(xv[1] - 1.0) < 1e-6) {
      bcv_model[v] = OPERATOR_BC_DIRICHLET;
      bcv_value[v] = exp(xv[0])*sin(xv[1]);
    }}

  Teuchos::RCP<PDE_SecondOrderPoisson> op_poisson = Teuchos::rcp(new PDE_SecondOrderPoisson(mesh));
  //The global stiffness operator
  Teuchos::RCP<Operator> global_op = op_poisson->global_operator();
  const CompositeVectorSpace & cvs = *op_poisson->GetCVS();
  //Before applying the boundary conditions we must populate the RHS
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("node");
  for (int v = 0; v < nnodes; v++) {
    mesh->node_get_coordinates(v, &xv);
    src[0][v] = 0;
  }
  global_op->UpdateRHS(source, true);
  op_poisson->SetBCs(bcv, bcv);
  //Global assembly
  op_poisson->ApplyBCs(true, true, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();
  //auto S = *global_op->A();
  //std::cout<<S<<std::endl;
  //Initialize the solution
  CompositeVector solution(cvs);
  solution.PutScalar(0.0);

  std::string xmlFileName = "test/operator_elasticity.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  Teuchos::ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
  global_op->InitializePreconditioner(slist);
  global_op->UpdatePreconditioner();

  Teuchos::ParameterList lop_list = plist.sublist("solvers")
                                         .sublist("PCG").sublist("pcg parameters");
  AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
      pcg(global_op, global_op);
  pcg.Init(lop_list);

  CompositeVector& rhs = *global_op->rhs();
  int ierr = pcg.ApplyInverse(rhs, solution);
  std::cout<<"sol is ="<<std::endl;
  solution.Print(std::cout);
  std::cout<<"rhs is"<<std::endl;
  rhs.Print(std::cout);
  std::cout<<"The solution is"<<std::endl;
  solution.Print(std::cout);
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh, true, false);
  io.InitializeCycle(0,0);
  const Epetra_MultiVector & sol = *solution.ViewComponent("node",false);
  std::cout<<"The epetra version of sol is"<<std::endl;
  std::cout<<sol<<std::endl;
  io.WriteVector(*sol(0), "solution", AmanziMesh::NODE);
  io.FinalizeCycle();
}