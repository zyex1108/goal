#include "mechanics.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "mesh.hpp"
#include "control.hpp"
#include "ev_bc_dirichlet.hpp"

namespace goal {

using Teuchos::rcpFromRef;

template <typename EvalT>
void goal::Mechanics::register_dirichlet(FieldManager fm)
{
  /* get the DBC parameterlist */
  RCP<const ParameterList> dbc_params =
    rcpFromRef(params->sublist("dirichlet bcs"));

  /* create a dummy data layout */
  RCP<Layouts> dl = rcp(new Layouts());

  /* temporary variable */
  RCP<PHX::Evaluator<GoalTraits> > ev;

  { /* dirichlet bc evaluator */
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<RCP<Layouts> >("Layouts", dl);
    p->set<RCP<Mesh> >("Mesh", mesh);
    p->set<RCP<Mechanics> >("Mechanics", rcpFromRef(*this));
    p->set<RCP<const ParameterList> >("DBC Parameters", dbc_params);
    ev = rcp(new BCDirichlet<EvalT, GoalTraits>(*p));
    fm->template registerEvaluator<EvalT>(ev);
    PHX::Tag<typename EvalT::ScalarT> tag("Dirichlet BCs", dl->dummy);
    fm->requireField<EvalT>(tag);
  }

}

}

template void goal::Mechanics::
register_dirichlet<goal::GoalTraits::Residual>(FieldManager fm);

template void goal::Mechanics::
register_dirichlet<goal::GoalTraits::Jacobian>(FieldManager fm);
