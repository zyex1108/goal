#include "ev_model_creep.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "workset.hpp"
#include "state_fields.hpp"
#include "control.hpp"
#include "assert_param.hpp"

namespace goal {

static RCP<ParameterList> get_valid_params()
{
  RCP<ParameterList> p = rcp(new ParameterList);
  p->set<double>("E", 0.0);
  p->set<double>("nu", 0.0);
  p->set<double>("K", 0.0);
  p->set<double>("Y", 0.0);
  p->set<double>("A", 0.0);
  p->set<double>("B", 0.0);
  p->set<double>("QR", 0.0);
  p->set<double>("C2", 0.0);
  return p;
}

static void validate_params(RCP<const ParameterList> p)
{
  assert_param(p, "E");
  assert_param(p, "nu");
  assert_param(p, "K");
  assert_param(p, "Y");
  assert_param(p, "A");
  assert_param(p, "B");
  assert_param(p, "QR");
  assert_param(p, "C2");
  p->validateParameters(*get_valid_params(), 0);
}

PHX_EVALUATOR_CTOR(ModelCreep, p) :
  dl            (p.get<RCP<Layouts> >("Layouts")),
  states        (p.get<RCP<StateFields> >("State Fields")),
  params        (p.get<RCP<const ParameterList> >("Material Params")),
  def_grad      (p.get<std::string>("Def Grad Name"), dl->qp_tensor),
  det_def_grad  (p.get<std::string>("Det Def Grad Name"), dl->qp_scalar),
  stress        (p.get<std::string>("Cauchy Name"), dl->qp_tensor)
{
  validate_params(params);
  E = params->get<double>("E");
  nu = params->get<double>("nu");
  K = params->get<double>("K");
  Y = params->get<double>("Y");
  A = params->get<double>("A");
  B = params->get<double>("B");
  QR = params->get<double>("QR");
  C2 = params->get<double>("C2");

  num_nodes = dl->node_qp_vector->dimension(1);
  num_qps = dl->node_qp_vector->dimension(2);
  num_dims = dl->node_qp_vector->dimension(3);

  this->addDependentField(def_grad);
  this->addDependentField(det_def_grad);
  this->addEvaluatedField(stress);
  this->setName("Model Creep");
}

PHX_POST_REGISTRATION_SETUP(ModelCreep, data, fm)
{
  this->utils.setFieldData(def_grad, fm);
  this->utils.setFieldData(det_def_grad, fm);
  this->utils.setFieldData(stress, fm);
}

PHX_EVALUATE_FIELDS(ModelCreep, workset)
{
}

}
