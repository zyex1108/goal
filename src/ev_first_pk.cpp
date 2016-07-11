#include "ev_first_pk.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "workset.hpp"
#include "phx_utils.hpp"

#include <Intrepid2_MiniTensor.h>

namespace goal {

PHX_EVALUATOR_CTOR(FirstPK, p) :
  dl            (p.get<RCP<Layouts> >("Layouts")),
  small_strain  (p.get<bool>("Small Strain")),
  def_grad      (p.get<std::string>("Def Grad Name"), dl->qp_tensor),
  det_def_grad  (p.get<std::string>("Det Def Grad Name"), dl->qp_scalar),
  cauchy        (p.get<std::string>("Cauchy Name"), dl->qp_tensor),
  first_pk      (p.get<std::string>("First PK Name"), dl->qp_tensor)
{
  num_qps = dl->node_qp_vector->dimension(2);
  num_dims = dl->node_qp_vector->dimension(3);

  this->addDependentField(def_grad);
  this->addDependentField(det_def_grad);
  this->addDependentField(cauchy);
  this->addEvaluatedField(first_pk);
  this->setName("First PK");
}

PHX_POST_REGISTRATION_SETUP(FirstPK, data, fm)
{
  this->utils.setFieldData(def_grad, fm);
  this->utils.setFieldData(det_def_grad, fm);
  this->utils.setFieldData(cauchy, fm);
  this->utils.setFieldData(first_pk, fm);
}

PHX_EVALUATE_FIELDS(FirstPK, workset)
{
  if (small_strain) {
    for (unsigned elem=0; elem < workset.size; ++elem)
    for (unsigned qp=0; qp < num_qps; ++qp)
    for (unsigned i=0; i < num_dims; ++i)
    for (unsigned j=0; j < num_dims; ++j)
      first_pk(elem,qp,i,j) = cauchy(elem,qp,i,j);
  }

  else {

    ScalarT J;
    Intrepid2::Tensor<ScalarT> F(num_dims);
    Intrepid2::Tensor<ScalarT> Finv(num_dims);
    Intrepid2::Tensor<ScalarT> sigma(num_dims);
    Intrepid2::Tensor<ScalarT> P(num_dims);
    Intrepid2::Tensor<ScalarT> I(Intrepid2::eye<ScalarT>(num_dims));

    for (unsigned elem=0; elem < workset.size; ++elem) {
      for (unsigned qp=0; qp < num_qps; ++qp) {

        J = det_def_grad(elem,qp);
        for (unsigned i=0; i < num_dims; ++i) {
          for (unsigned j=0; j < num_dims; ++j) {
            F(i,j) = def_grad(elem,qp,i,j);
            sigma(i,j) = cauchy(elem,qp,i,j);
          }
        }

        Finv = Intrepid2::inverse(F);
        ScalarT J = Intrepid2::det(F);
        P = J*sigma*Intrepid2::transpose(Finv);

        for (unsigned i=0; i < num_dims; ++i)
        for (unsigned j=0; j < num_dims; ++j)
          first_pk(elem,qp,i,j) = P(i,j);

      }
    }
  }
}

GOAL_INSTANTIATE_ALL(FirstPK)

}
