#include "solution_attachment.hpp"
#include "mesh.hpp"
#include "mechanics.hpp"
#include "solution_info.hpp"
#include "control.hpp"

#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>

namespace goal {

static Teuchos::Array<std::string> get_dual_names(RCP<Mechanics> mech)
{
  Teuchos::Array<std::string> names = mech->get_dof_names();
  Teuchos::Array<std::string> dual_names(0);
  for (unsigned i=0; i < names.size(); ++i)
    dual_names.push_back(names[i] + "_dual");
  return dual_names;
}

static void attach_vector_to_mesh(
    RCP<const Vector> u,
    RCP<Mesh> mesh,
    RCP<Mechanics> mech,
    Teuchos::Array<std::string> const& names,
    Teuchos::Array<std::string> const& offset_names)
{
  ArrayRCP<const ST> data = u->get1dView();
  apf::Mesh* m = mesh->get_apf_mesh();
  apf::DynamicArray<apf::Node> nodes = mesh->get_apf_nodes();
  std::vector<apf::Field*> fields;
  for (unsigned j=0; j < mesh->get_num_eqs(); ++j) {
    apf::Field* f = apf::createFieldOn(m, names[j].c_str(), apf::SCALAR);
    fields.push_back(f);
  }
  for (unsigned i=0; i < nodes.size(); ++i) {
    apf::Node* node = &(nodes[i]);
    if (! m->isOwned(node->entity)) continue;
    if (m->getType(node->entity) != apf::Mesh::VERTEX) continue;
    for (unsigned j=0; j < names.size(); ++j) {
      unsigned eq = mech->get_offset(offset_names[j]);
      LO row = mesh->get_lid(node, eq);
      double v = data[row];
      apf::setScalar(fields[j], node->entity, node->node, v);
    }
  }
  for (unsigned j=0; j < names.size(); ++j)
    apf::synchronize(fields[j]);
}

void attach_solutions_to_mesh(AttachInfo& ai)
{
  RCP<MultiVector> sv = ai.sol_info->owned_solution;
  unsigned nv = sv->getNumVectors();
  for (unsigned i=0; i < nv; ++i) {
    RCP<const Vector> u = sv->getVector(i);
    Teuchos::Array<std::string> names= ai.mech->get_var_names(i);
    attach_vector_to_mesh(u, ai.mesh, ai.mech, names, names);
  }
}

static void attach_vector_to_shape(
    RCP<const Vector> u,
    RCP<Mesh> mesh,
    RCP<Mechanics> mech,
    Teuchos::Array<std::string> const& names)
{
  ArrayRCP<const ST> data = u->get1dView();
  apf::Mesh* m = mesh->get_apf_mesh();
  apf::FieldShape* s = mesh->get_apf_shape();
  apf::DynamicArray<apf::Node> nodes = mesh->get_apf_nodes();
  std::vector<apf::Field*> fields;
  for (unsigned j=0; j < mesh->get_num_eqs(); ++j) {
    apf::Field* f = apf::createField(m, names[j].c_str(), apf::SCALAR, s);
    fields.push_back(f);
  }
  for (unsigned i=0; i < nodes.size(); ++i) {
    apf::Node* node = &(nodes[i]);
    if (! m->isOwned(node->entity)) continue;
    for (unsigned j=0; j < names.size(); ++j) {
      unsigned eq = mech->get_offset(names[j]);
      LO row = mesh->get_lid(node, eq);
      double v = data[row];
      apf::setScalar(fields[j], node->entity, node->node, v);
    }
  }
  for (unsigned j=0; j < names.size(); ++j)
    apf::synchronize(fields[j]);
}

void attach_solutions_to_shape(AttachInfo& ai)
{
  RCP<MultiVector> sv = ai.sol_info->owned_solution;
  unsigned nv = sv->getNumVectors();
  for (unsigned i=0; i < nv; ++i) {
    RCP<const Vector> u = sv->getVector(i);
    Teuchos::Array<std::string> names = ai.mech->get_var_names(i);
    attach_vector_to_shape(u, ai.mesh, ai.mech, names);
  }
}

void fill_vector_from_fields(
    RCP<Vector> u,
    RCP<Mesh> mesh,
    RCP<Mechanics> mech,
    Teuchos::Array<std::string> const& names)
{
  ArrayRCP<ST> data = u->get1dViewNonConst();
  apf::Mesh* m = mesh->get_apf_mesh();
  apf::DynamicArray<apf::Node> nodes = mesh->get_apf_nodes();
  std::vector<apf::Field*> fields;
  for (unsigned j=0; j < mesh->get_num_eqs(); ++j) {
    apf::Field* f = m->findField(names[j].c_str());
    CHECK(f);
    fields.push_back(f);
  }
  for (unsigned i=0; i < nodes.size(); ++i) {
    apf::Node* node = &(nodes[i]);
    if (! m->isOwned(node->entity)) continue;
    for (unsigned j=0; j < names.size(); ++j) {
      unsigned eq = mech->get_offset(names[j]);
      LO row = mesh->get_lid(node, eq);
      double v = apf::getScalar(fields[j], node->entity, node->node);
      data[row] = v;
    }
  }
}

void fill_solutions_from_fields(AttachInfo& ai)
{
  RCP<MultiVector> sv = ai.sol_info->owned_solution;
  unsigned nv = sv->getNumVectors();
  for (unsigned i=0; i < nv; ++i) {
    RCP<Vector> u = sv->getVectorNonConst(i);
    Teuchos::Array<std::string> names = ai.mech->get_var_names(i);
    fill_vector_from_fields(u, ai.mesh, ai.mech, names);
  }
}

void attach_dual_solutions_to_mesh(AttachInfo& ai)
{
  if (ai.sol_info->owned_dual == Teuchos::null) return;
  RCP<Vector> z = ai.sol_info->owned_dual;
  Teuchos::Array<std::string> offset_names = ai.mech->get_dof_names();
  Teuchos::Array<std::string> dual_names = get_dual_names(ai.mech);
  attach_vector_to_mesh(z, ai.mesh, ai.mech, dual_names, offset_names);
}

void remove_solutions_from_mesh(AttachInfo& ai)
{
  apf::Mesh* m = ai.mesh->get_apf_mesh();
  RCP<MultiVector> sv = ai.sol_info->owned_solution;
  unsigned nv = sv->getNumVectors();
  for (unsigned i=0; i < nv; ++i) {
    Teuchos::Array<std::string> names = ai.mech->get_var_names(i);
    for (unsigned j=0; j < names.size(); ++j) {
      apf::Field* f = m->findField(names[j].c_str());
      CHECK(f);
      apf::destroyField(f);
    }
  }
}

void remove_dual_solutions_from_mesh(AttachInfo& ai)
{
  if (ai.sol_info->owned_dual == Teuchos::null) return;
  apf::Mesh* m = ai.mesh->get_apf_mesh();
  Teuchos::Array<std::string> names = get_dual_names(ai.mech);
  for (unsigned i=0; i < names.size(); ++i) {
    apf::Field* f = m->findField(names[i].c_str());
    CHECK(f);
    apf::destroyField(f);
  }
}

}
