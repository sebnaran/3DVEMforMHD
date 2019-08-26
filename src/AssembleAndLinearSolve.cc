//Trilinos
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
//Amanzi
#include "VerboseObject_objs.hh"
#include "bilinear_form_registration.hh"
#include "MeshFactory.hh"
#include "LinearOperatorPCG.hh"
//#include "PDE_Elasticity.hh"
//Proper
#include "PDE_Test.hh"

int main(int argc, char *argv[]){
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;
  
  Teuchos::GlobalMPISession mpiSession(&argc, &argv); 
  auto comm = Amanzi::getDefaultComm();
  std::cout<<"got here1"<<std::endl;
  MeshFactory meshfactory(comm);
  //meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));

  // Generates a structured mesh covering [x0,x1] X [y0,y1] X [z0,z1] with
  // (nx, ny, nz) cells.
  // Teuchos::RCP<Mesh> create(const double x0, const double y0, const double z0,
  //                           const double x1, const double y1, const double z1,
  //                           const int nx, const int ny, const int nz,
  //                           const bool request_faces=true,
  //                           const bool request_edges=false);

  // Selected [-1,1] x[-1,1]x[-1,1] with nx =2 ny=2 nz =2
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0,2,2,2);
  std::cout<<"got here2"<<std::endl;
  //Teuchos::RCP<PDE_TestMHD> top = Teuchos::rcp(new PDE_TestMHD(mesh));


  // std::string xmlFileName = "test/operator_elasticity.xml";
  // Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  // Teuchos::ParameterList plist = xmlreader.getParameters();
  // Teuchos::ParameterList op_list = plist.sublist("PK operator")
  //                                       .sublist("elasticity operator");
  //std::cout<<plist<<std::endl;
  //Teuchos::RCP<PDE_Elasticity> op = Teuchos::rcp(new PDE_Elasticity(op_list, mesh));
  std::cout<<"got here3"<<std::endl;

  std::string xmlFileName = "test/operator_elasticity.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  Teuchos::ParameterList op_list = plist.sublist("PK operator")
                                        .sublist("elasticity operator");
  Teuchos::ParameterList& schema_list = op_list.sublist("schema");
  
  Teuchos::RCP<PDE_Test> op = Teuchos::rcp(new PDE_Test(mesh,op_list));
 
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  //int nnodes = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  //int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  //int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  //int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);
  



  Teuchos::RCP<Op> tLocal_Op = op->local_op();
  for( int c = 0 ; c < ncells ; c++){
    WhetStone::DenseMatrix M = tLocal_Op->matrices[c];
    std::cout<<"The matrix is"<<std::endl;
    PrintMatrix(M);
  }

  Teuchos::RCP<Operator> tGlobal_Op = op -> global_operator();
  //std::cout<<"got here4"<<std::endl;
  tGlobal_Op->SymbolicAssembleMatrix();
  std::cout<<"Got here5"<<std::endl;
  tGlobal_Op->AssembleMatrix();
  auto rcpA = tGlobal_Op->A();
  std::cout<<"The assembled matrix is"<<std::endl;
  std::cout<<*(rcpA.get())<<std::endl;
  CompositeVector& rhs = *tGlobal_Op ->rhs();
  auto & cvs = *(op->GetCVS());
  std::cout<<"Got here6"<<std::endl;
  std::cout<<"The rhs is"<<std::endl;
  rhs.Print(std::cout);

  Teuchos::ParameterList lop_list = plist.sublist("solvers")
                                         .sublist("PCG").sublist("pcg parameters");
  AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
        pcg(tGlobal_Op, tGlobal_Op);
  std::cout<<"Got here7"<<std::endl;
   CompositeVector solution(cvs);
   solution.PutScalar(0.0);
   pcg.Init(lop_list);
  std::cout<<"Got here8"<<std::endl;
  std::cout<<lop_list<<std::endl;
  std::cout<<"The initialized solution is"<<std::endl;
  solution.Print(std::cout);
    Teuchos::ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
  tGlobal_Op->InitializePreconditioner(slist);
  tGlobal_Op->UpdatePreconditioner();
  int ierr = pcg.ApplyInverse(rhs, solution);
  std::cout<<ierr<<std::endl;
  std::cout<<"The Solution is"<<std::endl;
  solution.Print(std::cout);
}