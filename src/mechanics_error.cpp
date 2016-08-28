#include "mechanics.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "mesh.hpp"
#include "control.hpp"
#include "assert_param.hpp"

#include "ev_gather_dual.hpp"
#include "ev_dof_interpolation.hpp"

static Teuchos::Array<std::string> get_dual_names(
    Teuchos::Array<std::string> const& names)
{
  Teuchos::Array<std::string> dual_names(0);
  for (unsigned i=0; i < names.size(); ++i)
    dual_names.push_back(names[i] + "_dual");
  return dual_names;
}

template <typename EvalT>
void goal::Mechanics::register_error(
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

  /* get the names of the dual variables */
  Teuchos::Array<std::string> names;
  get_dual_names(get_dof_names());

  { /* gather and interpolate the dual vector */
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<RCP<Layouts> >("Layouts", dl);
    p->set<RCP<Mesh> >("Mesh", mesh);
    p->set<Teuchos::Array<std::string> >("Dual Names", names);
    p->set<std::string>("BF Name", "BF");
    ev = rcp(new GatherDual<EvalT, GoalTraits>(*p));
    fm->template registerEvaluator<EvalT>(ev);
  }

  { /* evaluate the mechanics residual error contribution */
  }

}

/* ETI */
template void goal::Mechanics::
register_error<goal::GoalTraits::Forward>(
    std::string const& set, FieldManager fm);

template void goal::Mechanics::
register_error<goal::GoalTraits::Derivative>(
    std::string const& set, FieldManager fm);
