#include "ev_kinematics.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "workset.hpp"
#include "phx_utils.hpp"
#include <Intrepid2_MiniTensor.h>

namespace goal {

PHX_EVALUATOR_CTOR(Kinematics, p) :
  dl            (p.get<RCP<Layouts> >("Layouts")),
  disp_names    (p.get<Teuchos::Array<std::string> >("Disp Names")),
  def_grad      (p.get<std::string>("Def Grad Name"), dl->qp_tensor),
  det_def_grad  (p.get<std::string>("Det Def Grad Name"), dl->qp_scalar)
{
  num_qps = dl->node_qp_vector->dimension(2);
  num_dims = dl->node_qp_vector->dimension(3);

  grad_u.resize(num_dims);
  for (unsigned i=0; i < num_dims; ++i) {
    get_grad_field(disp_names[i], dl, grad_u[i]);
    this->addDependentField(grad_u[i]);
  }

  this->addEvaluatedField(def_grad);
  this->addEvaluatedField(det_def_grad);
  this->setName("Kinematics");
}

PHX_POST_REGISTRATION_SETUP(Kinematics, data, fm)
{
  for (unsigned i=0; i < num_dims; ++i)
    this->utils.setFieldData(grad_u[i], fm);
  this->utils.setFieldData(def_grad, fm);
  this->utils.setFieldData(det_def_grad, fm);
}

PHX_EVALUATE_FIELDS(Kinematics, workset)
{
  Intrepid2::Tensor<ScalarT> F(num_dims);
  for (unsigned elem=0; elem < workset.size; ++elem) {
    for (unsigned qp=0; qp < num_qps; ++qp) {
      for (unsigned i=0; i < num_dims; ++i)
      for (unsigned j=0; j < num_dims; ++j)
        def_grad(elem,qp,i,j) = grad_u[i](elem,qp,j);
      for (unsigned i=0; i < num_dims; ++i)
        def_grad(elem,qp,i,i) += 1.0;
      for (unsigned i=0; i < num_dims; ++i)
      for (unsigned j=0; j < num_dims; ++j)
        F(i,j) = def_grad(elem,qp,i,j);
      det_def_grad(elem,qp) = Intrepid2::det(F);
    }
  }
}

GOAL_INSTANTIATE_ALL(Kinematics)

}
