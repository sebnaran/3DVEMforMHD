
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

class PDE_FirstOrderPoisson: public PDE_HelperDiscretization {
 public:
//Constructor
  PDE_FirstOrderPoisson(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);//,Teuchos::ParameterList& plist);
//Destructor
~PDE_FirstOrderPoisson(){};
//   // -- creation of an operator
   using PDE_HelperDiscretization::UpdateMatrices;
   virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                               const Teuchos::Ptr<const CompositeVector>& p) override;

//   // -- postprocessing: calculated stress u from displacement p
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u) override {};

//These are accessors
Teuchos::RCP<CompositeVectorSpace> GetuCVS(){return u_cvs_;};
Teuchos::RCP<CompositeVectorSpace> GetpCVS(){return p_cvs_;};
//
Teuchos::RCP<Operator> GetMassGlobalOp(){return mass_global_op_;};
Teuchos::RCP<Operator> GetDivGlobalOp(){return div_global_op_;};
//
Teuchos::RCP<Op> GetMassLocalOp(){return mass_local_op_;};
Teuchos::RCP<Op> GetDivLocalOp(){return div_local_op_;};
//The Space of Velocities 
Teuchos::RCP<CompositeVectorSpace> u_cvs_;
//The Space of Pressures
Teuchos::RCP<CompositeVectorSpace> p_cvs_;
//Containers of the local matrices
Teuchos::RCP<Op> mass_local_op_;
Teuchos::RCP<Op> div_local_op_;
//Global Operators
Teuchos::RCP<Operator> mass_global_op_;
Teuchos::RCP<Operator> div_global_op_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif