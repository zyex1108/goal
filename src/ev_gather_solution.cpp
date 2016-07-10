#include "ev_gather_solution.hpp"
#include "mesh.hpp"
#include "workset.hpp"
#include "layouts.hpp"
#include "phx_utils.hpp"
#include "control.hpp"

namespace goal {

const std::string sol_names[3] =
{ "Solution",
  "Solution Dot",
  "Solution Dot Dot"};

template <typename Traits>
GatherSolution<GoalTraits::Residual, Traits>::
GatherSolution(ParameterList const& p) :
  dl      (p.get<RCP<Layouts> >("Layouts")),
  mesh    (p.get<RCP<Mesh> >("Mesh")),
  names   (p.get<Teuchos::Array<std::string> >("Sol Names")),
  index   (p.get<unsigned>("Sol Index"))
{
  num_eqs = names.size();
  num_nodes = dl->node_scalar->dimension(1);

  u.resize(num_eqs);
  for (unsigned i=0; i < num_eqs; ++i) {
    get_field<ScalarT>(names[i], dl, u[i]);
    this->addEvaluatedField(u[i]);
  }

  this->setName("Gather " + sol_names[index]);
}

template <typename Traits>
void GatherSolution<GoalTraits::Residual, Traits>::
postRegistrationSetup(
    typename Traits::SetupData d,
    PHX::FieldManager<Traits>& fm)
{
  for (unsigned i=0; i < num_eqs; ++i)
    this->utils.setFieldData(u[i], fm);
}

template <typename Traits>
void GatherSolution<GoalTraits::Residual, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  CHECK(workset.u != Teuchos::null);
  CHECK(workset.u->getNumVectors() >= index);
  RCP<const Vector> v = workset.u->getVector(index);
  ArrayRCP<const ST> sol = v->get1dView();
  CHECK(sol != Teuchos::null);

  for (unsigned elem=0; elem < workset.size; ++elem) {
    apf::MeshEntity* e = workset.ents[elem];
    for (unsigned node=0; node < num_nodes; ++node) {
      for (unsigned eq=0; eq < num_eqs; ++eq) {
        LO lid = mesh->get_lid(e, node, eq);
        u[eq](elem, node) = sol[lid];
      }
    }
  }
}

template <typename Traits>
GatherSolution<GoalTraits::Jacobian, Traits>::
GatherSolution(ParameterList const& p) :
  dl      (p.get<RCP<Layouts> >("Layouts")),
  mesh    (p.get<RCP<Mesh> >("Mesh")),
  names   (p.get<Teuchos::Array<std::string> >("Sol Names")),
  index   (p.get<unsigned>("Sol Index"))
{
  num_eqs = names.size();
  num_nodes = dl->node_scalar->dimension(1);

  u.resize(num_eqs);
  for (unsigned i=0; i < num_eqs; ++i) {
    get_field<ScalarT>(names[i], dl, u[i]);
    this->addEvaluatedField(u[i]);
  }

  this->setName("Gather " + sol_names[index]);
}

template <typename Traits>
void GatherSolution<GoalTraits::Jacobian, Traits>::
postRegistrationSetup(
    typename Traits::SetupData d,
    PHX::FieldManager<Traits>& fm)
{
  for (unsigned i=0; i < num_eqs; ++i)
    this->utils.setFieldData(u[i], fm);
}

template <typename Traits>
void GatherSolution<GoalTraits::Jacobian, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  CHECK(workset.u != Teuchos::null);
  CHECK(workset.u->getNumVectors() >= index);
  RCP<const Vector> v = workset.u->getVector(index);
  ArrayRCP<const ST> sol = v->get1dView();
  CHECK(sol != Teuchos::null);

  double fad_init = 0.0;
  if (index == 0) fad_init = workset.gamma;
  else if (index == 1) fad_init = workset.beta;
  else if (index == 2) fad_init = workset.alpha;

  for (unsigned elem=0; elem < workset.size; ++elem) {
    apf::MeshEntity* e = workset.ents[elem];
    for (unsigned node=0; node < num_nodes; ++node) {
      for (unsigned eq=0; eq < num_eqs; ++eq) {
        LO lid = mesh->get_lid(e, node, eq);
        unsigned offset = node*num_eqs + eq;
        u[eq](elem, node).val() = sol[lid];
        u[eq](elem, node).fastAccessDx(offset) = fad_init;
      }
    }
  }
}

GOAL_INSTANTIATE_ALL(GatherSolution)

}
