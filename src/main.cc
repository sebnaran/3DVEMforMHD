// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "exceptions.hh"
#include "Tensor.hh"
#include "CompositeVector.hh"

// Amanzi::Operators
#include "BilinearForm.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_HelperDiscretization.hh"
#include "Schema.hh"

#include <iostream>
#include <vector>

#include "UnitTest++.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

#include "MeshFactory.hh"
#include "Mesh_simple.hh"
#include "CompositeVector.hh"

#include "header.hh"

int main(){
    int a = 1, b = 2;
    int c;
    c = add(a,b);
    std::cout<<c<<std::endl;
    return 0;
}