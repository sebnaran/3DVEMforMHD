#include "Teuchos_GlobalMPISession.hpp"
#include "VerboseObject_objs.hh"
#include "bilinear_form_registration.hh"
#include "MeshFactory.hh"

int main(int argc, char *argv[]){
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv); 
  auto comm = Amanzi::getDefaultComm();
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 4, 5);
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  //  for (int c = 0; c < ncells; c++) {
  //    const Point& xc = mesh->cell_centroid(c);
  //  }
}