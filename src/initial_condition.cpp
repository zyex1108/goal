#include "initial_condition.hpp"
#include "mesh.hpp"
#include "mechanics.hpp"
#include "solution_info.hpp"
#include "expression.hpp"
#include "assert_param.hpp"
#include "control.hpp"

#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>

namespace goal {

static void validate_params(RCP<const ParameterList> p)
{
  assert_param(p, "ic");
  assert_param(p, "ic dot");
}

void fill_vector(
    RCP<Vector> vec,
    RCP<Mesh> mesh,
    Teuchos::Array<std::string> const& ic)
{
  ArrayRCP<ST> x = vec->get1dViewNonConst();
  apf::DynamicArray<apf::Node> nodes = mesh->get_apf_nodes();
  apf::Vector3 p;
  for (unsigned i=0; i < nodes.size(); ++i) {
    apf::Node node = nodes[i];
    mesh->get_apf_mesh()->getPoint(node.entity, node.node, p);
    for (unsigned eq=0; eq < mesh->get_num_eqs(); ++eq) {
      LO lid = mesh->get_lid(&node, eq);
      double v = expression_eval(ic[eq], p[0], p[1], p[2], 0.0);
      x[lid] = v;
    }
  }
}

void set_initial_conditions(
    RCP<const ParameterList> p,
    RCP<Mesh> mesh,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> info)
{
  double t0 = time();
  validate_params(p);
  Teuchos::Array<std::string> ic;
  Teuchos::Array<std::string> ic_dot;
  ic = p->get<Teuchos::Array<std::string> >("ic");
  ic_dot = p->get<Teuchos::Array<std::string> >("ic dot");
  CHECK(ic.size() == mech->get_num_eqs());
  CHECK(ic_dot.size() == mech->get_num_eqs());
  RCP<Vector> u = info->owned_solution->getVectorNonConst(0);
  RCP<Vector> v = info->owned_solution->getVectorNonConst(1);
  RCP<Vector> a = info->owned_solution->getVectorNonConst(2);
  fill_vector(u, mesh, ic);
  fill_vector(v, mesh, ic_dot);
  a->putScalar(0.0);
  info->scatter_solution();
  double t1 = time();
  print("initial conditions applied in %f seconds", t1-t0);
}

}
