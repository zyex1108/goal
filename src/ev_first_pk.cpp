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
  have_pressure (p.get<bool>("Have Pressure")),
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

  if (have_pressure) {
    std::string pressure_name = p.get<std::string>("Pressure Name");
    get_field(pressure_name, dl, pressure);
    this->addDependentField(pressure);
  }

  this->addEvaluatedField(first_pk);
  this->setName("First PK");
}

PHX_POST_REGISTRATION_SETUP(FirstPK, data, fm)
{
  this->utils.setFieldData(def_grad, fm);
  this->utils.setFieldData(det_def_grad, fm);
  this->utils.setFieldData(cauchy, fm);
  if (have_pressure)
    this->utils.setFieldData(pressure, fm);
  this->utils.setFieldData(first_pk, fm);
}

PHX_EVALUATE_FIELDS(FirstPK, workset)
{
  /* populate first pk tensor with cauchy tensor */
  for (unsigned elem=0; elem < workset.size; ++elem)
  for (unsigned qp=0; qp < num_qps; ++qp)
  for (unsigned i=0; i < num_dims; ++i)
  for (unsigned j=0; j < num_dims; ++j)
    first_pk(elem,qp,i,j) = cauchy(elem,qp,i,j);

  /* add in pressure if this is a mixed formulation */
  if (have_pressure) {
    for (unsigned elem=0; elem < workset.size; ++elem) {
      for (unsigned qp=0; qp < num_qps; ++qp) {
        ScalarT p = first_pk(elem, qp, 0, 0);
        for (unsigned i=1; i < num_dims; ++i)
          p += first_pk(elem, qp, i, i);
        p /= num_dims;
        for (unsigned i=0; i < num_dims; ++i)
          first_pk(elem,qp,i,i) += pressure(elem,qp) - p;
      }
    }
  }

  /* pull back to the reference config if this is large strain */
  if (! small_strain) {
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
