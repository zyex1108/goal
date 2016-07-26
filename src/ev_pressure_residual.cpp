#include "ev_pressure_residual.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "mesh.hpp"
#include "workset.hpp"
#include "phx_utils.hpp"
#include "assert_param.hpp"

#include <Intrepid2_MiniTensor.h>

namespace goal {

static void validate_params(RCP<const ParameterList> p)
{
  assert_param(p, "E");
  assert_param(p, "nu");
}


PHX_EVALUATOR_CTOR(PressureResidual, p) :
  dl            (p.get<RCP<Layouts> >("Layouts")),
  params        (p.get<RCP<const ParameterList> >("Material Params")),
  small_strain  (p.get<bool>("Small Strain")),
  pressure_name (p.get<std::string>("Pressure Name")),
  size          (p.get<std::string>("Size Name"), dl->qp_scalar),
  wDv           (p.get<std::string>("Weighted Dv Name"), dl->qp_scalar),
  BF            (p.get<std::string>("BF Name"), dl->node_qp_scalar),
  def_grad      (p.get<std::string>("Def Grad Name"), dl->qp_tensor),
  det_def_grad  (p.get<std::string>("Det Def Grad Name"), dl->qp_scalar),
  stress        (p.get<std::string>("Cauchy Name"), dl->qp_tensor)
{
  num_nodes = dl->node_qp_vector->dimension(1);
  num_qps = dl->node_qp_vector->dimension(2);
  num_dims = dl->node_qp_vector->dimension(3);

  validate_params(params);
  double E = params->get<double>("E");
  double nu = params->get<double>("nu");
  G = E/(2.0*(1.0+nu));
  K = E/(3.0*(1.0-2.0*nu));
  alpha = 1.0;

  get_grad_field(BF, dl, gBF);
  get_field(pressure_name, dl, pressure);
  get_grad_field(pressure_name, dl, pressure_grad);
  get_resid_field(pressure_name, dl, resid);

  this->addDependentField(size);
  this->addDependentField(wDv);
  this->addDependentField(BF);
  this->addDependentField(gBF);
  this->addDependentField(def_grad);
  this->addDependentField(det_def_grad);
  this->addDependentField(stress);
  this->addDependentField(pressure);
  this->addDependentField(pressure_grad);
  this->addEvaluatedField(resid);
  this->setName("Pressure Residual");
}

PHX_POST_REGISTRATION_SETUP(PressureResidual, data, fm)
{
  this->utils.setFieldData(size, fm);
  this->utils.setFieldData(wDv, fm);
  this->utils.setFieldData(BF, fm);
  this->utils.setFieldData(gBF, fm);
  this->utils.setFieldData(def_grad, fm);
  this->utils.setFieldData(det_def_grad, fm);
  this->utils.setFieldData(stress, fm);
  this->utils.setFieldData(pressure, fm);
  this->utils.setFieldData(pressure_grad, fm);
  this->utils.setFieldData(resid, fm);
}

PHX_EVALUATE_FIELDS(PressureResidual, workset)
{
  Intrepid2::Tensor<ScalarT> sigma(num_dims);
  Intrepid2::Tensor<ScalarT> F(num_dims);
  Intrepid2::Tensor<ScalarT> Cinv(num_dims);

  if (small_strain) {

    for (unsigned elem=0; elem < workset.size; ++elem) {
      for (unsigned node=0; node < num_nodes; ++node)
        resid(elem, node) = 0.0;
      for (unsigned qp=0; qp < num_qps; ++qp) {
        for (unsigned i=0; i < num_dims; ++i)
        for (unsigned j=0; j < num_dims; ++j)
          sigma(i,j) = stress(elem,qp,i,j);
        ScalarT dUdJ = (1.0/num_dims)*Intrepid2::trace(sigma);
        for (unsigned node=0; node < num_nodes; ++node)
          resid(elem,node) +=
            wDv(elem,qp)*BF(elem,node,qp)*(dUdJ - pressure(elem,qp))/K;
      }

      for (unsigned qp=0; qp < num_qps; ++qp) {
        double h = size(elem, qp);
        ScalarT param = 0.5*alpha*h*h / G;
        for (unsigned node=0; node < num_nodes; ++node)
        for (unsigned i=0; i < num_dims; ++i)
          resid(elem, node) -= param * wDv(elem,qp) *
            gBF(elem,node,qp,i)*pressure_grad(elem,qp,i);
      }
    }
  }

  else {

    for (unsigned elem=0; elem < workset.size; ++elem) {
      for (unsigned node=0; node < num_nodes; ++node)
        resid(elem, node) = 0.0;
      for (unsigned qp=0; qp < num_qps; ++qp) {
        for (unsigned i=0; i < num_dims; ++i)
        for (unsigned j=0; j < num_dims; ++j)
          sigma(i,j) = stress(elem,qp,i,j);
        ScalarT dUdJ = (1.0/num_dims)*Intrepid2::trace(sigma);
        for (unsigned node=0; node < num_nodes; ++node)
          resid(elem,node) +=
            wDv(elem,qp)*BF(elem,node,qp)*(dUdJ - pressure(elem,qp))/K;
      }

      for (unsigned qp=0; qp < num_qps; ++qp) {
        ScalarT J = det_def_grad(elem, qp);
        for (unsigned i=0; i < num_dims; ++i)
        for (unsigned j=0; j < num_dims; ++j)
          F(i,j) = def_grad(elem,qp,i,j);
        Cinv = Intrepid2::inverse(Intrepid2::transpose(F)*F);
        double h = size(elem, qp);
        ScalarT param = 0.5*alpha*h*h / G;
        for (unsigned node=0; node < num_nodes; ++node)
        for (unsigned i=0; i < num_dims; ++i)
        for (unsigned j=0; j < num_dims; ++j)
          resid(elem,node) -= param * wDv(elem, qp) * J * Cinv(i,j) *
            pressure_grad(elem,qp,i) * gBF(elem,node,qp,j);
      }
    }
  }

}

GOAL_INSTANTIATE_ALL(PressureResidual)

}
