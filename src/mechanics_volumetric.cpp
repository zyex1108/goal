#include "mechanics.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "mesh.hpp"
#include "ev_basis_functions.hpp"

template <typename EvalT>
void goal::Mechanics::register_volumetric(
    std::string const& set, FieldManager fm)
{

  /* do some work to create a data layout */
  unsigned ws_size = mesh->get_ws_size();
  unsigned n_nodes = mesh->get_num_elem_nodes();
  unsigned n_qps = mesh->get_num_elem_qps();
  unsigned n_dims = mesh->get_num_dims();
  RCP<Layouts> dl = rcp(new Layouts(ws_size, n_nodes, n_qps, n_dims));

  /* temporary variable */
  RCP<PHX::Evaluator<GoalTraits> > ev;

  { /* basis functions */
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<RCP<Layouts> >("Layouts", dl);
    p->set<RCP<Mesh> >("Mesh", mesh);
    p->set<std::string>("Weighted Dv Name", "wDv");
    p->set<std::string>("BF Name", "BF");
    ev = rcp(new BasisFunctions<EvalT, GoalTraits>(*p));
    fm->template registerEvaluator<EvalT>(ev);
  }

}

/* ETI */

template void goal::Mechanics::
register_volumetric<goal::GoalTraits::Residual>(
    std::string const& set, FieldManager fm);

template void goal::Mechanics::
register_volumetric<goal::GoalTraits::Jacobian>(
    std::string const& set, FieldManager fm);
