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
  /* parameters */
  ScalarT kappa = E/(3.0*(1.0-2.0*nu));
  ScalarT mu = E/(2.0*(1.0+nu));
  ScalarT sq23(std::sqrt(2.0/3.0));

  /* quantities at previous time */
  ScalarT eqps;
  Intrepid2::Tensor<ScalarT> Fp(num_dims);
  Intrepid2::Tensor<ScalarT> Fpinv(num_dims);
  Intrepid2::Tensor<ScalarT> Cpinv(num_dims);

  /* quantities at current time */
  ScalarT J;
  ScalarT Jm23;
  ScalarT dgam;
  ScalarT dgamp;
  Intrepid2::Tensor<ScalarT> F(num_dims);
  Intrepid2::Tensor<ScalarT> Fpn(num_dims);
  Intrepid2::Tensor<ScalarT> N(num_dims);
  Intrepid2::Tensor<ScalarT> sigma(num_dims);
  Intrepid2::Tensor<ScalarT> A(num_dims);
  Intrepid2::Tensor<ScalarT> expA(num_dims);
  Intrepid2::Tensor<ScalarT> I(Intrepid2::eye<ScalarT>(num_dims));

  /* trial state quantities */
  ScalarT f;
  ScalarT a0;
  ScalarT a1;
  ScalarT mubar;
  Intrepid2::Tensor<ScalarT> be(num_dims);
  Intrepid2::Tensor<ScalarT> s(num_dims);

  for (unsigned elem=0; elem < workset.size; ++elem) {

    apf::MeshEntity* e = workset.ents[elem];

    for (unsigned qp=0; qp < num_qps; ++qp) {

      /* compute the temperature adjusted relaxation parameter */
      B = A*std::exp(-QR / 303.0);

      /* deformation gradient quantities */
      for (unsigned i=0; i < num_dims; ++i)
      for (unsigned j=0; j < num_dims; ++j)
        F(i,j) = def_grad(elem, qp, i ,j);
      J = det_def_grad(elem, qp);
      Jm23 = std::pow(J, -2.0/3.0);

      /* get the plastic part of the def grad quantities */
      states->get_tensor("Fp_old", e, qp, Fp);
      Fpinv = Intrepid2::inverse(Fp);

      /* compute the trial state */
      Cpinv = Fpinv*Intrepid2::transpose(Fpinv);
      be = Jm23*F*Cpinv*Intrepid2::transpose(F);
      s = mu*Intrepid2::dev(be);
      mubar = Intrepid2::trace(be)*mu/num_dims;

      /* compute plasticity yield criteria */
      ScalarT smag = Intrepid2::norm<ScalarT>(s);
      states->get_scalar("eqps_old", e, qp, eqps);
      f = smag - sq23*(Y + K*eqps);

      /* compute creep onset criteria */
      a0 = Intrepid2::norm(Intrepid2::dev(be));
      a1 = Intrepid2::trace(be);
    }
  }


}

}
