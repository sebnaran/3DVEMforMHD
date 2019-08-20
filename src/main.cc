//Trilinos
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
//Amanzi
#include "VerboseObject_objs.hh"
#include "bilinear_form_registration.hh"
#include "MeshFactory.hh"
//#include "PDE_Elasticity.hh"
//Proper
#include "PDE_TestMHD.hh"

int main(int argc, char *argv[]){
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;
  
  Teuchos::GlobalMPISession mpiSession(&argc, &argv); 
  auto comm = Amanzi::getDefaultComm();

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
  
  //Teuchos::RCP<PDE_TestMHD> top = Teuchos::rcp(new PDE_TestMHD(mesh));


  std::string xmlFileName = "test/operator_elasticity.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  Teuchos::ParameterList op_list = plist.sublist("PK operator")
                                        .sublist("elasticity operator");
  std::cout<<plist<<std::endl;
  //Teuchos::RCP<PDE_Elasticity> op = Teuchos::rcp(new PDE_Elasticity(op_list, mesh));
  Teuchos::RCP<Test_MHD> op = Teuchos::rcp(new Test_MHD(mesh));
}