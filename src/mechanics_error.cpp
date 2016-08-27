#include "mechanics.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "mesh.hpp"
#include "control.hpp"
#include "assert_param.hpp"

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

}

/* ETI */
template void goal::Mechanics::
register_error<goal::GoalTraits::Forward>(
    std::string const& set, FieldManager fm);

template void goal::Mechanics::
register_error<goal::GoalTraits::Derivative>(
    std::string const& set, FieldManager fm);
