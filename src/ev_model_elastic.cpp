#include "ev_model_elastic.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "workset.hpp"
#include "state_fields.hpp"
#include "phx_utils.hpp"
#include "expression.hpp"
#include "assert_param.hpp"

namespace goal {

static RCP<ParameterList> get_valid_params()
{
  RCP<ParameterList> p = rcp(new ParameterList);
  p->set<double>("E", 0.0);
  p->set<double>("nu", 0.0);
  p->set<double>("alpha", 0.0);
  return p;
}

static void validate_params(
    RCP<const ParameterList> p,
    RCP<const ParameterList> tp)
{
  assert_param(p, "E");
  assert_param(p, "nu");
  if (tp != Teuchos::null)
    assert_param(p, "alpha");
  p->validateParameters(*get_valid_params(), 0);
}

PHX_EVALUATOR_CTOR(ModelElastic, p) :
  dl          (p.get<RCP<Layouts> >("Layouts")),
  states      (p.get<RCP<StateFields> >("State Fields")),
  params      (p.get<RCP<const ParameterList> >("Material Params")),
  temp_params (p.get<RCP<const ParameterList> >("Temperature Params")),
  disp_names  (p.get<Teuchos::Array<std::string> >("Disp Names")),
  stress      (p.get<std::string>("Cauchy Name"), dl->qp_tensor)
{
  validate_params(params, temp_params);
  E = params->get<double>("E");
  nu = params->get<double>("nu");

  have_temp = (Teuchos::nonnull(temp_params));
  if (have_temp) alpha = params->get<double>("alpha");

  num_nodes = dl->node_qp_vector->dimension(1);
  num_qps = dl->node_qp_vector->dimension(2);
  num_dims = dl->node_qp_vector->dimension(3);

  grad_u.resize(num_dims);
  for (unsigned i=0; i < num_dims; ++i) {
    get_grad_field(disp_names[i], dl, grad_u[i]);
    this->addDependentField(grad_u[i]);
  }

  this->addEvaluatedField(stress);
  this->setName("Model Elastic");
}

PHX_POST_REGISTRATION_SETUP(ModelElastic, data, fm)
{
  for (unsigned i=0; i < num_dims; ++i)
    this->utils.setFieldData(grad_u[i], fm);
  this->utils.setFieldData(stress, fm);
}

PHX_EVALUATE_FIELDS(ModelElastic, workset)
{
  ScalarT lambda = E*nu/((1.0+nu)*(1.0-2.0*nu));
  ScalarT mu = E/(2.0*(1.0+nu));

  Intrepid2::Tensor<ScalarT> eps(num_dims);
  Intrepid2::Tensor<ScalarT> sigma(num_dims);
  Intrepid2::Tensor<ScalarT> I(Intrepid2::eye<ScalarT>(num_dims));

  for (unsigned elem=0; elem < workset.size; ++elem) {
    for (unsigned qp=0; qp < num_qps; ++qp) {

      for (unsigned i=0; i < num_dims; ++i)
      for (unsigned j=0; j < num_dims; ++j)
        eps(i,j) = 0.5*(grad_u[i](elem,qp,j) + grad_u[j](elem,qp,i));

      sigma = eps*(2.0*mu) + I*(lambda*Intrepid2::trace(eps));

      for (unsigned i=0; i < num_dims; ++i)
      for (unsigned j=0; j < num_dims; ++j)
        stress(elem,qp,i,j) = sigma(i,j);

      apf::MeshEntity* e = workset.ents[elem];
      states->set_tensor("cauchy",e,qp,sigma);

    }
  }

  if (have_temp) {

    double three_kappa = E/(1.0-2.0*nu);
    std::string val = temp_params->get<std::string>("value");
    double T = expression_eval(val,0,0,0,workset.t_new);
    double T_ref = temp_params->get<double>("reference");

    for (unsigned elem=0; elem < workset.size; ++elem)
    for (unsigned qp=0; qp < num_qps; ++qp)
    for (unsigned i=0; i < num_dims; ++i)
      stress(elem,qp,i,i) -= three_kappa*alpha*(T-T_ref);
  }

}

GOAL_INSTANTIATE_ALL(ModelElastic)

}
