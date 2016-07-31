#include "ev_qoi_avg_displacement.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "workset.hpp"
#include "phx_utils.hpp"

namespace goal {

PHX_EVALUATOR_CTOR(QoIAvgDisplacement, p) :
  dl          (p.get<RCP<Layouts> >("Layouts")),
  disp_names  (p.get<Teuchos::Array<std::string> >("Disp Names")),
  wDv         (p.get<std::string>("Weighted Dv Name"), dl->qp_scalar),
  avg_disp    (p.get<std::string>("Avg Displacement Name"), dl->elem_scalar)
{
  num_nodes = dl->node_qp_vector->dimension(1);
  num_qps = dl->node_qp_vector->dimension(2);
  num_dims = dl->node_qp_vector->dimension(3);

  disp.resize(disp_names.size());
  for (unsigned i=0; i < disp_names.size(); ++i) {
    get_field(disp_names[i], dl, disp[i]);
    this->addDependentField(disp[i]);
  }

  this->addDependentField(wDv);
  this->addEvaluatedField(avg_disp);
  this->setName("QoI: Average Displacement");
}

PHX_POST_REGISTRATION_SETUP(QoIAvgDisplacement, data, fm)
{
  this->utils.setFieldData(wDv, fm);
  for (unsigned i=0; i < disp_names.size(); ++i)
    this->utils.setFieldData(disp[i], fm);
  this->utils.setFieldData(avg_disp, fm);
}

PHX_EVALUATE_FIELDS(QoIAvgDisplacement, workset)
{
  for (unsigned elem=0; elem < workset.size; ++elem) {
    avg_disp(elem) = 0.0;
    for (unsigned qp=0; qp < num_qps; ++qp)
    for (unsigned i=0; i < num_dims; ++i)
      avg_disp(elem) += disp[i](elem, qp) * wDv(elem, qp);
    avg_disp(elem) /= num_dims;
  }
}

GOAL_INSTANTIATE_JACOBIAN(QoIAvgDisplacement)

}
