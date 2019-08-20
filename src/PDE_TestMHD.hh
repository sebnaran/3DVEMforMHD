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

class Test_MHD: public PDE_HelperDiscretization {
 public:

  Test_MHD(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh):
  PDE_HelperDiscretization(mesh)
  {
    global_op_ = Teuchos::null;
    //operator_type_ = OPERATOR_ELASTICITY;
  }

//   PDE_Elasticity(Teuchos::ParameterList& plist,
//                  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
//       PDE_HelperDiscretization(mesh),
//       K_(Teuchos::null),
//       K_default_(1.0)
//   {
//     global_op_ = Teuchos::null;
//     operator_type_ = OPERATOR_ELASTICITY;
//     Init_(plist);
//   }

//   // main virtual members
//   // -- setup 
//   void SetTensorCoefficient(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K);
//   void SetTensorCoefficient(double K);

//   // -- creation of an operator
   using PDE_HelperDiscretization::UpdateMatrices;
   virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                               const Teuchos::Ptr<const CompositeVector>& p) override;

//   // -- postprocessing: calculated stress u from displacement p
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u) override {};

//   // access
//   const Schema& global_schema_col() { return global_schema_col_; }
//   const Schema& schema_col() { return local_schema_col_; }
//   const Schema& schema_row() { return local_schema_row_; }

//  protected:
//   void Init_(Teuchos::ParameterList& plist);

//  protected:
//   Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;
//   double K_default_;

//   Teuchos::RCP<WhetStone::BilinearForm> mfd_;

//   // operator and schemas
//   Schema global_schema_col_, global_schema_row_;
//   Schema local_schema_col_, local_schema_row_;

//   OperatorType operator_type_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif