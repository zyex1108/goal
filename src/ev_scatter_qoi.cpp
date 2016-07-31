#include "ev_scatter_qoi.hpp"
#include "mesh.hpp"
#include "control.hpp"
#include "workset.hpp"
#include "layouts.hpp"

namespace goal {

template <typename Traits>
ScatterQoI<GoalTraits::Residual, Traits>::
ScatterQoI(ParameterList const& p) :
  dl      (p.get<RCP<Layouts> >("Layouts")),
  qoi     (p.get<std::string>("QoI Name"), dl->elem_scalar)
{
  this->addDependentField(qoi);
  std::string name = "Scatter QoI";
  PHX::Tag<ScalarT> op(name, dl->dummy);
  this->setName("name");
}

template <typename Traits>
void ScatterQoI<GoalTraits::Residual, Traits>::
postRegistrationSetup(
    typename Traits::SetupData d,
    PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(qoi, fm);
}

template <typename Traits>
void ScatterQoI<GoalTraits::Residual, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
}

template <typename Traits>
ScatterQoI<GoalTraits::Jacobian, Traits>::
ScatterQoI(ParameterList const& p) :
  dl      (p.get<RCP<Layouts> >("Layouts")),
  mesh    (p.get<RCP<Mesh> >("Mesh")),
  qoi     (p.get<std::string>("QoI Name"), dl->elem_scalar)
{
  this->addDependentField(qoi);
  std::string name = "Scatter QoI";
  PHX::Tag<ScalarT> op(name, dl->dummy);
  this->setName("name");
}

template <typename Traits>
void ScatterQoI<GoalTraits::Jacobian, Traits>::
postRegistrationSetup(
    typename Traits::SetupData d,
    PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(qoi, fm);
}

template <typename Traits>
void ScatterQoI<GoalTraits::Jacobian, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
}

GOAL_INSTANTIATE_ALL(ScatterQoI)

}
