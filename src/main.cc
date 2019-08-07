#include "UnitTest++.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

#include "MeshFactory.hh"
#include "Mesh_simple.hh"
#include "CompositeVector.hh"
#include "header.hh"

using namespace Teuchos;
int main(){
    int a = 1, b = 2;
    int c;
    c = add(a,b);
    std::cout<<c<<std::endl;
    return 0;
}