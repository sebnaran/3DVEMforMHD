
#ifndef AMANZI_OPERATOR_PDE_TESTMHD_HH_
#define AMANZI_OPERATOR_PDE_TESTMHD_HH_

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

namespace Amanzi {
namespace Operators {

class PDE_ElectroMagMHD: public PDE_HelperDiscretization {
 public:
//Constructor
  PDE_ElectroMagMHD(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)://,Teuchos::ParameterList& plist);
  PDE_HelperDiscretization(mesh){}
//Destructor
~PDE_ElectroMagMHD(){};
//   // -- creation of an operator
   using PDE_HelperDiscretization::UpdateMatrices;
   virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                               const Teuchos::Ptr<const CompositeVector>& p) override;

//   // -- postprocessing: calculated stress u from displacement p
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u) override {};


void LocalMagMatrix();
void ContructLocalMatrix();
//containers for local matrices
Teuchos::RCP<Op> mag_local_op_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif