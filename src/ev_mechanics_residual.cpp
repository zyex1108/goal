#include "ev_mechanics_residual.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "workset.hpp"
#include "phx_utils.hpp"

namespace goal {

PHX_EVALUATOR_CTOR(MechanicsResidual, p) :
  dl          (p.get<RCP<Layouts> >("Layouts")),
  bf_name     (p.get<std::string>("BF Name")),
  disp_names  (p.get<Teuchos::Array<std::string> >("Disp Names")),
  wDv         (p.get<std::string>("Weighted Dv Name"), dl->qp_scalar),
  stress      (p.get<std::string>("First PK Name"), dl->qp_tensor)
{
  num_nodes = dl->node_qp_vector->dimension(1);
  num_qps = dl->node_qp_vector->dimension(2);
  num_dims = dl->node_qp_vector->dimension(3);

  resid.resize(num_dims);
  for (unsigned i=0; i < num_dims; ++i) {
    get_resid_field(disp_names[i], dl, resid[i]);
    this->addEvaluatedField(resid[i]);
  }

  get_grad_field(bf_name, dl, gBF);

  this->addDependentField(wDv);
  this->addDependentField(gBF);
  this->addDependentField(stress);
  this->setName("Mechanics Residual");
}

PHX_POST_REGISTRATION_SETUP(MechanicsResidual, data, fm)
{
  this->utils.setFieldData(wDv, fm);
  this->utils.setFieldData(gBF, fm);
  this->utils.setFieldData(stress, fm);
  for (unsigned i=0; i < num_dims; ++i)
    this->utils.setFieldData(resid[i], fm);
}

PHX_EVALUATE_FIELDS(MechanicsResidual, workset)
{
  for (unsigned elem=0; elem < workset.size; ++elem) {
    for (unsigned node=0; node < num_nodes; ++node)
    for (unsigned dim=0; dim < num_dims; ++dim)
      resid[dim](elem, node) = ScalarT(0.0);
    for (unsigned qp=0; qp < num_qps; ++qp)
    for (unsigned node=0; node < num_nodes; ++node)
    for (unsigned i=0; i < num_dims; ++i)
    for (unsigned j=0; j < num_dims; ++j)
      resid[i](elem, node) +=
        stress(elem,qp,i,j)*gBF(elem,node,qp,j)*wDv(elem,qp);
  }
}

GOAL_INSTANTIATE_ALL(MechanicsResidual)

}
