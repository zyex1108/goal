#include "ev_bc_dirichlet.hpp"
#include "layouts.hpp"
#include "workset.hpp"
#include "mesh.hpp"
#include "mechanics.hpp"
#include "expression.hpp"
#include "control.hpp"

#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>

namespace goal {

template <typename Traits>
void BCDirichlet<GoalTraits::Residual, Traits>::
validate_params()
{
  using Teuchos::Array;
  using Teuchos::ParameterList;
  using Teuchos::ParameterEntry;
  using Teuchos::getValue;

  for (auto it=params->begin(); it != params->end(); ++it) {
    ParameterEntry const& entry = params->entry(it);
    Array<std::string> a = getValue<Array<std::string> >(entry);
    CHECK(a.size() == 3);
    std::string const& dof = a[0];
    std::string const& set = a[1];
    mesh->get_nodes(set);
    mechanics->get_offset(dof);
  }
}

static double get_bc_val(
    std::string const& val,
    RCP<Mesh> m,
    apf::Node* n,
    const double t)
{
  double v = 0.0;
  apf::Mesh* apfm = m->get_apf_mesh();
  apf::MeshEntity* ent = n->entity;
  if (apfm->getType(ent) == apf::Mesh::VERTEX) {
    apf::Vector3 p;
    apfm->getPoint(ent, 0, p);
    v = expression_eval(val, p[0], p[1], p[2], t);
  }
  return v;
}

template <typename Traits>
BCDirichlet<GoalTraits::Residual, Traits>::
BCDirichlet(ParameterList const& p) :
  dl        (p.get<RCP<Layouts> >("Layouts")),
  mesh      (p.get<RCP<Mesh> >("Mesh")),
  mechanics (p.get<RCP<Mechanics> >("Mechanics")),
  params    (p.get<RCP<const ParameterList> >("DBC Parameters"))
{
  validate_params();

  std::string name = "Dirichlet BCs";
  PHX::Tag<ScalarT> op(name, dl->dummy);
  this->setName(name);
  this->addEvaluatedField(op);
}

template <typename Traits>
void BCDirichlet<GoalTraits::Residual, Traits>::
postRegistrationSetup(
    typename Traits::SetupData d,
    PHX::FieldManager<Traits>& fm)
{
}

template <typename Traits>
void BCDirichlet<GoalTraits::Residual, Traits>::
apply_bc(
    typename Traits::EvalData workset,
    Teuchos::Array<std::string> const& a)
{
  CHECK(workset.r != Teuchos::null);
  CHECK(workset.u != Teuchos::null);

  std::string const& dof = a[0];
  std::string const& set = a[1];
  std::string const& val = a[2];

  apf::Vector3 p;
  unsigned offset = mechanics->get_offset(dof);

  ArrayRCP<ST> res = workset.r->get1dViewNonConst();
  RCP<const Vector> u = workset.u->getVector(0);
  ArrayRCP<const ST> sol = u->get1dView();
  CHECK(sol != Teuchos::null);
  CHECK(res != Teuchos::null);

  double t = workset.t_new;
  std::vector<apf::Node*> const& nodes = mesh->get_nodes(set);
  for (unsigned i=0; i < nodes.size(); ++i) {
    apf::Node* node = nodes[i];
    LO row = mesh->get_lid(node, offset);
    double v = get_bc_val(val, mesh, node, t);
    res[row] = sol[row] - v;
  }
}

template <typename Traits>
void BCDirichlet<GoalTraits::Residual, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  using Teuchos::Array;
  using Teuchos::ParameterList;
  using Teuchos::ParameterEntry;
  using Teuchos::getValue;

  for (auto i=this->params->begin(); i != this->params->end(); ++i) {
    ParameterEntry const& entry = this->params->entry(i);
    Array<std::string> a = getValue<Array<std::string> >(entry);
    this->template apply_bc(workset, a);
  }
}

template <typename Traits>
void BCDirichlet<GoalTraits::Jacobian, Traits>::
validate_params()
{
  using Teuchos::Array;
  using Teuchos::ParameterList;
  using Teuchos::ParameterEntry;
  using Teuchos::getValue;

  for (auto it=params->begin(); it != params->end(); ++it) {
    ParameterEntry const& entry = params->entry(it);
    Array<std::string> a = getValue<Array<std::string> >(entry);
    CHECK(a.size() == 3);
    std::string const& dof = a[0];
    std::string const& set = a[1];
    mesh->get_nodes(set);
    mechanics->get_offset(dof);
  }
}

template <typename Traits>
BCDirichlet<GoalTraits::Jacobian, Traits>::
BCDirichlet(ParameterList const& p) :
  dl        (p.get<RCP<Layouts> >("Layouts")),
  mesh      (p.get<RCP<Mesh> >("Mesh")),
  mechanics (p.get<RCP<Mechanics> >("Mechanics")),
  params    (p.get<RCP<const ParameterList> >("DBC Parameters"))
{
  validate_params();

  std::string name = "Dirichlet BCs";
  PHX::Tag<ScalarT> op(name, dl->dummy);
  this->setName(name);
  this->addEvaluatedField(op);
}

template <typename Traits>
void BCDirichlet<GoalTraits::Jacobian, Traits>::
postRegistrationSetup(
    typename Traits::SetupData d,
    PHX::FieldManager<Traits>& fm)
{
}

template <typename Traits>
void BCDirichlet<GoalTraits::Jacobian, Traits>::
apply_bc(
    typename Traits::EvalData workset,
    Teuchos::Array<std::string> const& a)
{
  CHECK(workset.J != Teuchos::null);

  std::string const& dof = a[0];
  std::string const& set = a[1];
  std::string const& val = a[2];

  bool fill_res = (workset.r != Teuchos::null);

  ArrayRCP<const ST> sol;
  ArrayRCP<ST> res;

  if (fill_res) {
    CHECK(workset.u != Teuchos::null);
    CHECK(workset.r != Teuchos::null);
    RCP<const Vector> u = workset.u->getVector(0);
    sol = u->get1dView();
    res = workset.r->get1dViewNonConst();
    CHECK(sol != Teuchos::null);
    CHECK(res != Teuchos::null);
  }

  Teuchos::Array<LO> index(1);
  Teuchos::Array<ST> value(1);
  value[0] = 1.0;
  size_t num_entries;
  Teuchos::Array<LO> matrix_indices;
  Teuchos::Array<ST> matrix_entries;

  unsigned offset = mechanics->get_offset(dof);
  RCP<Matrix> J = workset.J;
  double t = workset.t_new;
  std::vector<apf::Node*> const& nodes = mesh->get_nodes(set);

  for (unsigned i=0; i < nodes.size(); ++i) {

    apf::Node* node = nodes[i];
    LO row = mesh->get_lid(node, offset);

    if (fill_res) {
      double v = get_bc_val(val, mesh, node, t);
      res[row] = sol[row] - v;
    }

    index[0] = row;
    num_entries = J->getNumEntriesInLocalRow(row);
    matrix_indices.resize(num_entries);
    matrix_entries.resize(num_entries);
    J->getLocalRowCopy(row, matrix_indices(), matrix_entries(), num_entries);
    for (unsigned c=0; c < num_entries; ++c)
      matrix_entries[c] = 0.0;
    J->replaceLocalValues(row, matrix_indices(), matrix_entries());
    J->replaceLocalValues(row, index(), value());
  }
}

template <typename Traits>
void BCDirichlet<GoalTraits::Jacobian, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  using Teuchos::Array;
  using Teuchos::ParameterList;
  using Teuchos::ParameterEntry;
  using Teuchos::getValue;

  for (auto i=this->params->begin(); i != this->params->end(); ++i) {
    ParameterEntry const& entry = this->params->entry(i);
    Array<std::string> a = getValue<Array<std::string> >(entry);
    this->template apply_bc(workset, a);
  }
}

GOAL_INSTANTIATE_ALL(BCDirichlet)

}
