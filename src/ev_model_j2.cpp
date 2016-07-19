#include "ev_model_j2.hpp"
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
  return p;
}

static void validate_params(RCP<const ParameterList> p)
{
  assert_param(p, "E");
  assert_param(p, "nu");
  assert_param(p, "K");
  assert_param(p, "Y");
  p->validateParameters(*get_valid_params(), 0);
}

PHX_EVALUATOR_CTOR(ModelJ2, p) :
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

  num_nodes = dl->node_qp_vector->dimension(1);
  num_qps = dl->node_qp_vector->dimension(2);
  num_dims = dl->node_qp_vector->dimension(3);

  this->addDependentField(def_grad);
  this->addDependentField(det_def_grad);
  this->addEvaluatedField(stress);
  this->setName("Model J2");
}

PHX_POST_REGISTRATION_SETUP(ModelJ2, data, fm)
{
  this->utils.setFieldData(def_grad, fm);
  this->utils.setFieldData(det_def_grad, fm);
  this->utils.setFieldData(stress, fm);
}

PHX_EVALUATE_FIELDS(ModelJ2, workset)
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
  Intrepid2::Tensor<ScalarT> F(num_dims);
  Intrepid2::Tensor<ScalarT> Fpn(num_dims);
  Intrepid2::Tensor<ScalarT> N(num_dims);
  Intrepid2::Tensor<ScalarT> sigma(num_dims);
  Intrepid2::Tensor<ScalarT> I(Intrepid2::eye<ScalarT>(num_dims));

  /* trial state quantities */
  ScalarT f;
  ScalarT mubar;
  Intrepid2::Tensor<ScalarT> be(num_dims);
  Intrepid2::Tensor<ScalarT> s(num_dims);

  for (unsigned elem=0; elem < workset.size; ++elem) {

    apf::MeshEntity* e = workset.ents[elem];

    for (unsigned qp=0; qp < num_qps; ++qp) {

      /* deformation gradient quantities */
      for (unsigned i=0; i < num_dims; ++i)
      for (unsigned j=0; j < num_dims; ++j)
        F(i,j) = def_grad(elem, qp, i, j);
      J = det_def_grad(elem, qp);
      Jm23 = std::pow(J, -2.0/3.0);

      /* get plastic part of def grad quantities */
      states->get_tensor("Fp_old", e, qp, Fp);
      Fpinv = Intrepid2::inverse(Fp);

      /* compute the trial state */
      Cpinv = Fpinv*Intrepid2::transpose(Fpinv);
      be = Jm23*F*Cpinv*Intrepid2::transpose(F);
      s = mu*Intrepid2::dev(be);
      mubar = Intrepid2::trace(be)*mu/num_dims;

      /* check yield condition */
      ScalarT smag = Intrepid2::norm<ScalarT>(s);
      states->get_scalar("eqps_old", e, qp, eqps);
      ScalarT f = smag - sq23 * (Y + K*eqps);

      /* plastic increment - return mapping algorithm */
      if (f > 1.0e-12) {

        bool converged = false;
        dgam = 0.0;
        ScalarT H = 0.0;
        ScalarT dH = 0.0;
        ScalarT alpha = 0.0;
        ScalarT res = 0.0;
        unsigned iter = 0;

        ScalarT X = 0.0;
        ScalarT R = f;
        ScalarT dRdX = -2.0*mubar*(1.0+H/(3.0*mubar));

        while (!converged && iter < 30) {
          iter++;
          X = X - R/dRdX;
          alpha = eqps + sq23*X;
          H = K*alpha;
          dH = K;
          R = smag - (2.0*mubar*X + sq23*(Y+H));
          dRdX = -2.0*mubar*(1.0+dH/(3.0*mubar));
          res = std::abs(R);
          if ((res < 1.0e-11) || (res/Y < 1.0e-11) || (res/f < 1.0e-11))
            converged = true;
          if (iter == 30)
            fail("J2: return mapping failed to converge");
        }

        /* updates */
        dgam = X;
        N = (1.0/smag)*s;
        s -= 2.0*mubar*dgam*N;
        states->set_scalar("eqps", e, qp, alpha);

        /* get Fpn */
        Fpn = Intrepid2::exp(dgam*N)*Fp;
        states->set_tensor("Fp", e, qp, Fpn);
      }

      /* otherwise elastic increment */
      else
        states->set_scalar("eqps", e, qp, eqps);

      /* compute stress */
      ScalarT p = 0.5*kappa*(J-1.0/J);
      sigma = I*p + s/J;
      states->set_tensor("cauchy", e, qp, sigma);
      for (unsigned i=0; i < num_dims; ++i)
      for (unsigned j=0; j < num_dims; ++j)
        stress(elem, qp, i, j) = sigma(i, j);
    }
  }
}

GOAL_INSTANTIATE_ALL(ModelJ2)

}
