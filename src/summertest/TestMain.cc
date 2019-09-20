#include <UnitTest++.h>

#include "elasticitytest.cc"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_GlobalMPISession.hpp"

#include "VerboseObject_objs.hh"
#include "bilinear_form_registration.hh"


int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return UnitTest::RunAllTests();
}
