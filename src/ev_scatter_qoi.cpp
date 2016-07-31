#include "ev_scatter_qoi.hpp"
#include "mesh.hpp"
#include "control.hpp"
#include "workset.hpp"
#include "layouts.hpp"

namespace goal {

using Teuchos::ArrayRCP;

template <typename Traits>
ScatterQoI<GoalTraits::Forward, Traits>::
ScatterQoI(ParameterList const& p) :
  dl      (p.get<RCP<Layouts> >("Layouts")),
  qoi     (p.get<std::string>("QoI Name"), dl->elem_scalar)
{
  std::string name = "Scatter QoI";
  PHX::Tag<ScalarT> op(name, dl->dummy);
  this->addDependentField(qoi);
  this->addEvaluatedField(op);
  this->setName("name");
}

template <typename Traits>
void ScatterQoI<GoalTraits::Forward, Traits>::
postRegistrationSetup(
    typename Traits::SetupData d,
    PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(qoi, fm);
}

template <typename Traits>
void ScatterQoI<GoalTraits::Forward, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
}

template <typename Traits>
ScatterQoI<GoalTraits::Derivative, Traits>::
ScatterQoI(ParameterList const& p) :
  dl      (p.get<RCP<Layouts> >("Layouts")),
  mesh    (p.get<RCP<Mesh> >("Mesh")),
  qoi     (p.get<std::string>("QoI Name"), dl->elem_scalar)
{


  std::string name = "Scatter QoI";
  PHX::Tag<ScalarT> op(name, dl->dummy);
  this->addDependentField(qoi);
  this->addEvaluatedField(op);
  this->setName("name");
}

template <typename Traits>
void ScatterQoI<GoalTraits::Derivative, Traits>::
postRegistrationSetup(
    typename Traits::SetupData d,
    PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(qoi, fm);
}

template <typename Traits>
void ScatterQoI<GoalTraits::Derivative, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  CHECK(workset.q != Teuchos::null);
  ArrayRCP<ST> dqdu = workset.q->get1dViewNonConst();
  CHECK(dqdu != Teuchos::null);

  unsigned num_eqs = mesh->get_num_eqs();
  for (unsigned elem=0; elem < workset.size; ++elem) {
    apf::MeshEntity* e = workset.ents[elem];
    unsigned dof = 0;
    for (unsigned node=0; node < num_nodes; ++node) {
      for (unsigned eq=0; eq < num_eqs; ++eq) {
        LO row = mesh->get_lid(e, node, eq);
        dqdu[row] += qoi(elem).fastAccessDx(dof);
        ++dof;
      }
    }
  }
}

GOAL_INSTANTIATE_ALL(ScatterQoI)

}
