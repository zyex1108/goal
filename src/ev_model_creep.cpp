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
  A2 = params->get<double>("A");
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

  /* time increment quantities */
  double dt = workset.t_new - workset.t_old;

  for (unsigned elem=0; elem < workset.size; ++elem) {

    apf::MeshEntity* e = workset.ents[elem];

    for (unsigned qp=0; qp < num_qps; ++qp) {

      /* compute the temperature adjusted relaxation parameter */
      B = A2*std::exp(-QR / 303.0);

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

      /* below yield strength */
      if (f <= 0.0) {

        /* creep increment - return mapping algorithm */
        if (a0 > 1.0e-12) {

          bool converged = false;
          ScalarT res = 0.0;
          unsigned iter = 0;

          ScalarT X = 1.1e-4;
          ScalarT R = X - dt*B*std::pow(mu, C2)*
            std::pow((a0 - 2.0/3.0*X*a1)*(a0 - 2.0/3.0*X*a1), C2/2.0);

          std::cout << R << std::endl;

          ScalarT dRdX = 1.0 - dt*B*std::pow(mu, C2)*(C2/2.0)*
            std::pow((a0 - 2.0/3.0*X*a1)*(a0 - 2.0/3.0*X*a1), C2/2.0 - 1.0)*
            (8.0/9.0*X*a1*a1 - 4.0/3.0*a0*a1);

          while (!converged && iter < 30) {
            iter++;
            X = X - R/dRdX;
            R = X - dt*B*std::pow(mu, C2)*
              std::pow((a0 - 2.0/3.0*X*a1)*(a0 - 2.0/3.0*X*a1), C2/2.0);
            dRdX = 1.0 - dt*B*std::pow(mu, C2)*(C2/2.0)*
              std::pow((a0 - 2.0/3.0*X*a1)*(a0 - 2.0/3.0*X*a1), C2/2.0 - 1.0)*
              (8.0/9.0*X*a1*a1 - 4.0/3.0*a0*a1);
            res = std::abs(R);
            if (res < 1.0e-10)
              converged = true;
            if (iter == 30)
              fail("Creep: pure creep increment failed to converge");
          }

          /* updates */
          dgam = X;
          N = (1.0/smag)*s;
          s -= 2.0*mubar*dgam*N;

          /* exponential map to get Fpnew */
          A = dgam*N;
          expA = Intrepid2::exp(A);
          Fpn = expA*Fp;
          states->set_tensor("Fp", e, qp, Fpn);
          states->set_scalar("eqps", e, qp, eqps);
        }

        /* purely elastic increment - no creep */
        else {
          states->set_scalar("eqps", e, qp, eqps);
        }
      }

      /* plastic increment - return mapping algorithm */
      else {

        bool converged = false;
        ScalarT H = 0.0;
        ScalarT dH = 0.0;
        ScalarT alpha = 0.0;
        ScalarT res = 0.0;
        unsigned iter = 0;

        dgam = 0.0;
        dgamp = 0.0;

        ScalarT X = 0.0;
        ScalarT R = f;
        ScalarT dRdX = -2.0*mubar*(1.0+H/(2.0*mubar));

        while (! converged) {
          iter++;
          X = X - R/dRdX;
          H = 2.0*mubar*dt*B*
            std::pow((smag+(2.0/3.0)*K*X-f)*(smag+(2.0/3.0)*K*X-f), C2/2.0);
          dH = (4.0/3.0)*C2*mubar*dt*B*K*
            std::pow((smag+(2.0/3.0)*K*X-f)*(smag+(2.0/3.0)*K*X-f), (C2-1.0)/2.0);
          R = f - 2.0*mubar*(1.0+K/(3.0*mubar))*X - H;
          dRdX = -2.0*mubar*(1.0+K/(3.0*mubar)) - dH;
          res = std::abs(R);
          if ((res < 1.0e-10) || (res/f) < 1.0e-11)
            converged = true;
          if (iter == 30)
            fail("Creep: plastic increment failed to converge");
        }

        /* updates */
        dgamp = X;
        N = s / Intrepid2::norm(s);
        s -= -2.0*mubar*dgamp*N + f*N - 2.0*mubar*(1.0+K/(3.0*mubar))*dgamp*N;
        dgam = dgamp + dt*B*std::pow(Intrepid2::norm(s), C2);
        alpha = eqps + sq23*dgamp;
        N = s / Intrepid2::norm(s);
        states->set_scalar("eqps", e, qp, alpha);

        /* exponential map to get Fpnew */
        A = dgam*N;
        expA = Intrepid2::exp(A);
        Fpn = expA*Fp;
        states->set_tensor("Fp", e, qp, Fpn);

      }

      /* compute stress */
      ScalarT p = 0.5*kappa*(J-1.0/J);
      sigma = I*p + s/J;
      states->set_tensor("cauchy", e, qp, sigma);
      for (unsigned i=0; i < num_dims; ++i)
      for (unsigned j=0; j < num_dims; ++j)
        stress(elem, qp, i, j) = sigma(i,j);

    }
  }
}

GOAL_INSTANTIATE_ALL(ModelCreep)

}
