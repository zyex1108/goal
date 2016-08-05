#include "ev_bc_traction.hpp"
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

template <typename EvalT, typename Traits>
void BCTraction<EvalT, Traits>::validate_params()
{
  using Teuchos::Array;
  using Teuchos::ParameterEntry;
  using Teuchos::getValue;

  for (auto it=params->begin(); it != params->end(); ++it) {
    ParameterEntry const& entry = params->entry(it);
    CHECK(entry.isType<Array<std::string> >());
    Array<std::string> a = getValue<Array<std::string> >(entry);
    CHECK(a.size() == mesh->get_num_dims() + 1);
    std::string const& set = a[0];
    mesh->get_facets(set);
  }
}

PHX_EVALUATOR_CTOR(BCTraction, p) :
  dl        (p.get<RCP<Layouts> >("Layouts")),
  mesh      (p.get<RCP<Mesh> >("Mesh")),
  mechanics (p.get<RCP<Mechanics> >("Mechanics")),
  params    (p.get<RCP<const ParameterList> >("NBC Parameters"))
{
  validate_params();
  apf_mesh = mesh->get_apf_mesh();
  num_dims = mesh->get_num_dims();

  std::string name = "Neumann BCs";
  PHX::Tag<ScalarT> op(name, dl->dummy);
  this->setName(name);
  this->addEvaluatedField(op);
}

PHX_POST_REGISTRATION_SETUP(BCTraction, data, fm)
{
}

template <typename EvalT, typename Traits>
void BCTraction<EvalT, Traits>::apply_bc(
    typename Traits::EvalData workset,
    Teuchos::Array<std::string> const& a)
{
  CHECK(workset.r != Teuchos::null);
  Teuchos::ArrayRCP<ST> res;
  res = workset.r->get1dViewNonConst();
  CHECK(res != Teuchos::null);

  std::string const& set = a[0];
  std::vector<apf::MeshEntity*> const& facets = mesh->get_facets(set);
  apf::GlobalNumbering* numbering = mesh->get_apf_numbering();
  unsigned q_order = mesh->get_q_order();
  unsigned num_qps = mesh->get_num_elem_qps();
  
  for (unsigned i=0; i < facets.size(); ++i) {
    apf::MeshEntity* f = facets[i];
    apf::MeshElement* me = apf::createMeshElement(apf_mesh, f);
    unsigned num_nodes = apf::getElementNumbers(numbering, f, numbers);
    for (unsigned qp=0; qp < num_qps; ++qp) {
      apf::getIntPoint(me, q_order, qp, local);
      apf::mapLocalToGlobal(me, local, global);
      apf::getBF(mesh->get_apf_shape(), me, local, BF);
      double w = apf::getIntWeight(me, q_order, qp);
      double dv = apf::getDV(me, local);
      for (unsigned i=0; i < num_dims; ++i)
        traction[i] = expression_eval(
            a[i+1], global[0], global[1], global[2], workset.t_new);
      for (unsigned node=0; node < num_nodes; ++node) {
        for (unsigned eq=0; eq < num_dims; ++eq) {
          LO lid = mesh->get_lid(f, node, eq);
          res[lid] -= BF[node]*traction[eq]*dv*w;
        }
      }
    }
  }
}

PHX_EVALUATE_FIELDS(BCTraction, workset)
{
  using Teuchos::Array;
  using Teuchos::ParameterEntry;
  using Teuchos::getValue;

  for (auto i=params->begin(); i != params->end(); ++i) {
    ParameterEntry const& entry = params->entry(i);
    Array<std::string> a = getValue<Array<std::string> >(entry);
    apply_bc(workset, a);
  }
}

GOAL_INSTANTIATE_ALL(BCTraction)

}
