#include "ev_pressure_residual.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "workset.hpp"
#include "phx_utils.hpp"

namespace goal {

PHX_EVALUATOR_CTOR(PressureResidual, p) :
  dl            (p.get<RCP<Layouts> >("Layouts")),
  pressure_name (p.get<std::string>("Pressure Name")),
  wDv           (p.get<std::string>("Weighted Dv Name"), dl->qp_scalar),
  BF            (p.get<std::string>("BF Name"), dl->node_qp_scalar),
  def_grad      (p.get<std::string>("Def Grad Name"), dl->qp_tensor),
  stress        (p.get<std::string>("First PK Name"), dl->qp_tensor)
{
  num_nodes = dl->node_qp_vector->dimension(1);
  num_qps = dl->node_qp_vector->dimension(2);
  num_dims = dl->node_qp_vector->dimension(3);

  get_grad_field(BF, dl, gBF);
  get_field(pressure_name, dl, pressure);
  get_grad_field(pressure_name, dl, pressure_grad);
  get_resid_field(pressure_name, dl, resid);

  this->addDependentField(wDv);
  this->addDependentField(BF);
  this->addDependentField(gBF);
  this->addDependentField(stress);
  this->addDependentField(pressure);
  this->addDependentField(pressure_grad);
  this->addEvaluatedField(resid);
  this->setName("Pressure Residual");
}

PHX_POST_REGISTRATION_SETUP(PressureResidual, data, fm)
{
  this->utils.setFieldData(wDv, fm);
  this->utils.setFieldData(BF, fm);
  this->utils.setFieldData(gBF, fm);
  this->utils.setFieldData(stress, fm);
  this->utils.setFieldData(pressure, fm);
  this->utils.setFieldData(pressure_grad, fm);
  this->utils.setFieldData(resid, fm);
}

PHX_EVALUATE_FIELDS(PressureResidual, workset)
{
}

GOAL_INSTANTIATE_ALL(PressureResidual)

}
