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

  Teuchos::RCP<Operator> Global_Op00 = op -> global_operator();
  //std::cout<<"got here4"<<std::endl;
  Global_Op00->SymbolicAssembleMatrix();
  std::cout<<"Got here5"<<std::endl;
  Global_Op00->AssembleMatrix();
  auto rcpA00 = Global_Op00->A();
  std::cout<<"The assembled matrix is"<<std::endl;
  std::cout<<*(rcpA00.get())<<std::endl;
  CompositeVector& rhs00 = *Global_Op00 ->rhs();
  const CompositeVectorSpace & cvs00 = *(op->GetCVS());
  std::cout<<"Got here6"<<std::endl;
  std::cout<<"The rhs is"<<std::endl;
  rhs00.Print(std::cout);
  
  Teuchos::RCP<Operator> Global_Op11 = op -> global_operator();
  Global_Op11->SymbolicAssembleMatrix();
  std::cout<<"Got here5"<<std::endl;
  Global_Op11->AssembleMatrix();
  auto rcpA11 = Global_Op11->A();
  std::cout<<"The assembled matrix is"<<std::endl;
  std::cout<<*(rcpA11.get())<<std::endl;
  CompositeVector& rhs11 = *Global_Op11 ->rhs();
  const auto & cvs11 = *(op->GetCVS());
  std::cout<<"Got here6"<<std::endl;
  std::cout<<"The rhs is"<<std::endl;
  rhs11.Print(std::cout);


  Teuchos::RCP<TreeVectorSpace> tvs = Teuchos::rcp(new TreeVectorSpace());
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(cvs00))));
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(cvs11))));

  Teuchos::RCP<TreeOperator> Treeop = Teuchos::rcp(new Operators::TreeOperator(tvs));
  Treeop->SetOperatorBlock(0,0,Global_Op00);
  Treeop->SetOperatorBlock(1,1,Global_Op11); 

  Treeop->SymbolicAssembleMatrix();
  Treeop->AssembleMatrix();


  std::cout<<"about to print the tree matrix"<<std::endl;
  auto A = *Treeop->A();
  std::cout<<"The tree matrix has been assembled"<<std::endl;
  std::cout<<A<<std::endl;
  
  TreeVector rhs(*tvs);
  *rhs.SubVector(0)->Data() = *Global_Op00->rhs();
  *rhs.SubVector(1)->Data() = *Global_Op11->rhs();
  std::cout<<"The right-hand-side is"<<std::endl;
  rhs.Print(std::cout);

  TreeVector solution(*tvs);
  solution.PutScalar(0.0);
   Teuchos::ParameterList lop_list = plist.sublist("solvers")
                                         .sublist("PCG").sublist("pcg parameters");
  AmanziSolvers::LinearOperatorPCG<TreeOperator, TreeVector, TreeVectorSpace>
                pcg(Treeop, Treeop);
  // std::cout<<"Got here7"<<std::endl;

    pcg.Init(lop_list);
    std::cout<<"Got here8"<<std::endl;
   std::cout<<lop_list<<std::endl;
   Teuchos::ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
   Global_Op00->InitializePreconditioner(slist);
   Global_Op00->UpdatePreconditioner();
   Global_Op11->InitializePreconditioner(slist);
   Global_Op11->UpdatePreconditioner();
   Treeop->InitBlockDiagonalPreconditioner();
   //Treeop->InitializePreconditioner(slist);
   //Treeop->UpdatePreconditioner();
   int ierr = pcg.ApplyInverse(rhs, solution);
   std::cout<<ierr<<std::endl;
   std::cout<<"The Solution is"<<std::endl;
   solution.Print(std::cout);
}