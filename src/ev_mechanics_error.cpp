#include "ev_mechanics_error.hpp"
#include "mesh.hpp"
#include "workset.hpp"
#include "layouts.hpp"
#include "phx_utils.hpp"
#include "control.hpp"

namespace goal {

PHX_EVALUATOR_CTOR(MechanicsError, p) :
  dl      (p.get<RCP<Layouts> >("Layouts")),
  mesh    (p.get<RCP<Mesh> >("Mesh"))
{
  num_nodes = dl->node_qp_scalar->dimension(1);
  num_qps = dl->node_qp_vector->dimension(2);
  num_dims = dl->node_qp_vector->dimension(3);

  this->setName("Mechanics Error");
}

PHX_POST_REGISTRATION_SETUP(MechanicsError, data, fm)
{
}

PHX_EVALUATE_FIELDS(MechanicsError, workset)
{
}

GOAL_INSTANTIATE_ALL(MechanicsError)

}
