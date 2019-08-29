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

//Proper
#include "PDE_FirstOrderPoisson.hh"

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
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(-1,-1,1,1,3,3);
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  Teuchos::RCP<PDE_FirstOrderPoisson> op_poisson = Teuchos::rcp(new PDE_FirstOrderPoisson(mesh));
  
  //The global mass matrix operator
  Teuchos::RCP<Operator> mass_op = op_poisson->GetMassGlobalOp();
  //Global assembly
  mass_op->SymbolicAssembleMatrix();
  mass_op->AssembleMatrix();
  auto M = mass_op->A();
  std::cout<<*M<<std::endl;

  Teuchos::RCP<Op> mass_local_op = op_poisson->GetMassLocalOp();
   for( int c = 0 ; c < ncells ; c++){
    WhetStone::DenseMatrix Mi = mass_local_op->matrices[c];
     std::cout<<"The matrix is"<<std::endl;
     std::cout<<Mi<<std::endl;
     //PrintMatrix(Mi);
   }
  //The global divergence matrix operator
  Teuchos::RCP<Operator> div_op = op_poisson->GetDivGlobalOp();
  //Global assembly
  div_op->SymbolicAssembleMatrix();
  div_op->AssembleMatrix();

  auto S = mass_op->A();
  std::cout<<*S<<std::endl;

  
  // Teuchos::RCP<PDE_Test> op = Teuchos::rcp(new PDE_Test(mesh,op_list));
 
  // int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  // //int nnodes = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  // //int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  // //int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  // //int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);
  



  // Teuchos::RCP<Op> tLocal_Op = op->local_op();
  // for( int c = 0 ; c < ncells ; c++){
  //   WhetStone::DenseMatrix M = tLocal_Op->matrices[c];
  //   std::cout<<"The matrix is"<<std::endl;
  //   PrintMatrix(M);
  // }

  // Teuchos::RCP<Operator> Global_Op00 = op -> global_operator();
  // //std::cout<<"got here4"<<std::endl;
  // Global_Op00->SymbolicAssembleMatrix();
  // std::cout<<"Got here5"<<std::endl;
  // Global_Op00->AssembleMatrix();
  // auto rcpA00 = Global_Op00->A();
  // std::cout<<"The assembled matrix is"<<std::endl;
  // std::cout<<*(rcpA00.get())<<std::endl;
  // CompositeVector& rhs00 = *Global_Op00 ->rhs();
  // const CompositeVectorSpace & cvs00 = *(op->GetCVS());
  // std::cout<<"Got here6"<<std::endl;
  // std::cout<<"The rhs is"<<std::endl;
  // rhs00.Print(std::cout);
  
  // Teuchos::RCP<Operator> Global_Op11 = op -> global_operator();
  // Global_Op11->SymbolicAssembleMatrix();
  // std::cout<<"Got here5"<<std::endl;
  // Global_Op11->AssembleMatrix();
  // auto rcpA11 = Global_Op11->A();
  // std::cout<<"The assembled matrix is"<<std::endl;
  // std::cout<<*(rcpA11.get())<<std::endl;
  // CompositeVector& rhs11 = *Global_Op11 ->rhs();
  // const auto & cvs11 = *(op->GetCVS());
  // std::cout<<"Got here6"<<std::endl;
  // std::cout<<"The rhs is"<<std::endl;
  // rhs11.Print(std::cout);


  // Teuchos::RCP<TreeVectorSpace> tvs = Teuchos::rcp(new TreeVectorSpace());
  // tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(cvs00))));
  // tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(cvs11))));

  // Teuchos::RCP<TreeOperator> Treeop = Teuchos::rcp(new Operators::TreeOperator(tvs));
  // Treeop->SetOperatorBlock(0,0,Global_Op00);
  // Treeop->SetOperatorBlock(1,1,Global_Op11); 

  // Treeop->SymbolicAssembleMatrix();
  // Treeop->AssembleMatrix();


  // std::cout<<"about to print the tree matrix"<<std::endl;
  // auto A = *Treeop->A();
  // std::cout<<"The tree matrix has been assembled"<<std::endl;
  // std::cout<<A<<std::endl;
  
  // TreeVector rhs(*tvs);
  // *rhs.SubVector(0)->Data() = *Global_Op00->rhs();
  // *rhs.SubVector(1)->Data() = *Global_Op11->rhs();
  // std::cout<<"The right-hand-side is"<<std::endl;
  // rhs.Print(std::cout);

  // TreeVector solution(*tvs);
  // solution.PutScalar(0.0);
  //  Teuchos::ParameterList lop_list = plist.sublist("solvers")
  //                                        .sublist("PCG").sublist("pcg parameters");
  // AmanziSolvers::LinearOperatorPCG<TreeOperator, TreeVector, TreeVectorSpace>
  //               pcg(Treeop, Treeop);
  // // std::cout<<"Got here7"<<std::endl;

  //   pcg.Init(lop_list);
  //   std::cout<<"Got here8"<<std::endl;
  //  std::cout<<lop_list<<std::endl;
  //  Teuchos::ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
  //  Global_Op00->InitializePreconditioner(slist);
  //  Global_Op00->UpdatePreconditioner();
  //  Global_Op11->InitializePreconditioner(slist);
  //  Global_Op11->UpdatePreconditioner();
  //  Treeop->InitBlockDiagonalPreconditioner();
  //  //Treeop->InitializePreconditioner(slist);
  //  //Treeop->UpdatePreconditioner();
  //  int ierr = pcg.ApplyInverse(rhs, solution);
  //  std::cout<<ierr<<std::endl;
  //  std::cout<<"The Solution is"<<std::endl;
  //  solution.Print(std::cout);
}