#include "mechanics.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "mesh.hpp"
#include "control.hpp"
#include "assert_param.hpp"

#include "ev_qoi_avg_displacement.hpp"
#include "ev_scatter_qoi.hpp"

template <typename EvalT>
void goal::Mechanics::register_qoi(
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

  /* get the quantity of interest parameters */
  assert_sublist(params, "qoi");
  RCP<const ParameterList> qoi_params = rcpFromRef(params->sublist("qoi"));
  std::string qoi_name = qoi_params->get<std::string>("name");

  if (qoi_name == "avg displacement") {
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<RCP<Layouts> >("Layouts", dl);
    p->set<Teuchos::Array<std::string> >("Disp Names", fields["disp"]);
    p->set<std::string>("Weighted Dv Name", "wDv");
    p->set<std::string>("Avg Displacement Name", qoi_name);
    ev = rcp(new QoIAvgDisplacement<EvalT, GoalTraits>(*p));
    fm->template registerEvaluator<EvalT>(ev);
  }

  else
    fail("unknown qoi name: %s", qoi_name.c_str());

  { /* scatter qoi */
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<RCP<Layouts> >("Layouts", dl);
    p->set<RCP<Mesh> >("Mesh", mesh);
    p->set<std::string>("QoI Name", qoi_name);
    ev = rcp(new ScatterQoI<EvalT, GoalTraits>(*p));
    fm->template registerEvaluator<EvalT>(ev);
    PHX::Tag<typename EvalT::ScalarT> tag("Scatter QoI", dl->dummy);
    fm->requireField<EvalT>(tag);
  }

}

/* ETI */
template void goal::Mechanics::
register_qoi<goal::GoalTraits::Forward>(
    std::string const& set, FieldManager fm);

template void goal::Mechanics::
register_qoi<goal::GoalTraits::Derivative>(
    std::string const& set, FieldManager fm);
