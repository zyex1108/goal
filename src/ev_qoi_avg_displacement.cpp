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
}

PHX_POST_REGISTRATION_SETUP(QoIAvgDisplacement, data, fm)
{
}

PHX_EVALUATE_FIELDS(QoIAvgDisplacement, workset)
{
}

GOAL_INSTANTIATE_JACOBIAN(QoIAvgDisplacement)

}
