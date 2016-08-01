#include "ev_scatter_residual.hpp"
#include "mesh.hpp"
#include "control.hpp"
#include "workset.hpp"
#include "layouts.hpp"
#include "phx_utils.hpp"

namespace goal {

using Teuchos::ArrayRCP;
using Teuchos::arrayView;

template <typename Traits>
ScatterResidual<GoalTraits::Forward, Traits>::
ScatterResidual(ParameterList const& p) :
  dl        (p.get<RCP<Layouts> >("Layouts")),
  mesh      (p.get<RCP<Mesh> >("Mesh")),
  dof_names (p.get<Teuchos::Array<std::string> >("DOF Names"))
{
  num_nodes = dl->node_vector->dimension(1);
  num_eqs = dof_names.size();

  resid.resize(num_eqs);
  for (unsigned i=0; i < num_eqs; ++i) {
    get_resid_field(dof_names[i], dl, resid[i]);
    this->addDependentField(resid[i]);
  }

  std::string name = "Scatter Residual";
  PHX::Tag<ScalarT> op(name, dl->dummy);
  this->addEvaluatedField(op);
  this->setName(name);
}

template <typename Traits>
void ScatterResidual<GoalTraits::Forward, Traits>::
postRegistrationSetup(
    typename Traits::SetupData d,
    PHX::FieldManager<Traits>& fm)
{
  for (unsigned i=0; i < num_eqs; ++i)
    this->utils.setFieldData(resid[i], fm);
}

template <typename Traits>
void ScatterResidual<GoalTraits::Forward, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  CHECK(workset.r != Teuchos::null);
  ArrayRCP<ST> r = workset.r->get1dViewNonConst();
  CHECK(r != Teuchos::null);

  for (unsigned elem=0; elem < workset.size; ++elem) {
    apf::MeshEntity* e= workset.ents[elem];
    for (unsigned node=0; node < num_nodes; ++node) {
      for (unsigned eq=0; eq < num_eqs; ++eq) {
        LO lid = mesh->get_lid(e, node, eq);
        r[lid] += resid[eq](elem, node);
      }
    }
  }
}

template <typename Traits>
ScatterResidual<GoalTraits::Derivative, Traits>::
ScatterResidual(ParameterList const& p) :
  dl        (p.get<RCP<Layouts> >("Layouts")),
  mesh      (p.get<RCP<Mesh> >("Mesh")),
  dof_names (p.get<Teuchos::Array<std::string> >("DOF Names"))
{
  num_nodes = dl->node_vector->dimension(1);
  num_eqs = dof_names.size();
  num_dofs = num_nodes * num_eqs;

  resid.resize(num_eqs);
  for (unsigned i=0; i < num_eqs; ++i) {
    get_resid_field(dof_names[i], dl, resid[i]);
    this->addDependentField(resid[i]);
  }

  std::string name = "Scatter Residual";
  PHX::Tag<ScalarT> op(name, dl->dummy);
  this->addEvaluatedField(op);
  this->setName(name);
}

template <typename Traits>
void ScatterResidual<GoalTraits::Derivative, Traits>::
postRegistrationSetup(
    typename Traits::SetupData d,
    PHX::FieldManager<Traits>& fm)
{
  for (unsigned i=0; i < num_eqs; ++i)
    this->utils.setFieldData(resid[i], fm);
}

template <typename Traits>
void ScatterResidual<GoalTraits::Derivative, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  CHECK(workset.J != Teuchos::null);
  RCP<Matrix> J = workset.J;

  bool fill_resid = false;
  ArrayRCP<ST> r;
  if (workset.r != Teuchos::null) {
    fill_resid = true;
    r = workset.r->get1dViewNonConst();
    CHECK(r != Teuchos::null);
  }

  Teuchos::Array<LO> cols(num_dofs);

  if (! workset.is_adjoint) {
    for (unsigned elem=0; elem < workset.size; ++elem) {
      apf::MeshEntity* e = workset.ents[elem];
      unsigned idx=0;
      for (unsigned node=0; node < num_nodes; ++node) {
        for (unsigned eq=0; eq < num_eqs; ++eq) {
          cols[idx] = mesh->get_lid(e, node, eq);
          idx++;
        }
      }
      for (unsigned node=0; node < num_nodes; ++node) {
        for (unsigned eq=0; eq < num_eqs; ++eq) {
          LO row = mesh->get_lid(e, node, eq);
          FadType v = resid[eq](elem, node);
          J->sumIntoLocalValues(
              row, cols, arrayView(&(v.fastAccessDx(0)), num_dofs));
          if (fill_resid)
            r[row] += resid[eq](elem, node).val();
        }
      }
    }
  }

  else {
    for (unsigned elem=0; elem < workset.size; ++elem) {
      apf::MeshEntity* e = workset.ents[elem];
      unsigned idx=0;
      for (unsigned node=0; node < num_nodes; ++node) {
        for (unsigned eq=0; eq < num_eqs; ++eq) {
          cols[idx] = mesh->get_lid(e, node, eq);
          idx++;
        }
      }
      for (unsigned node=0; node < num_nodes; ++node) {
        for (unsigned eq=0; eq < num_eqs; ++eq) {
          LO row = mesh->get_lid(e, node, eq);
          FadType v = resid[eq](elem, node);
          for (unsigned dof=0; dof < num_dofs; ++dof)
            J->sumIntoLocalValues(cols[dof], arrayView(&row, 1),
                arrayView(&(v.fastAccessDx(dof)), 1));
        }
      }
    }
  }
}

GOAL_INSTANTIATE_ALL(ScatterResidual)

}
