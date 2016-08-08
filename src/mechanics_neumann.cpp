#include "mechanics.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "mesh.hpp"
#include "ev_bc_traction.hpp"

namespace goal {

using Teuchos::rcpFromRef;

template <typename EvalT>
void goal::Mechanics::register_neumann(FieldManager fm)
{

  /* do nothing if there are no neumann bcs */
  if (! params->isSublist("neumann bcs")) return;

  /* validate parameters */
  RCP<const ParameterList> nbc_params =
    rcpFromRef(params->sublist("neumann bcs"));
  RCP<ParameterList> vp = rcp(new ParameterList);
  vp->sublist("tractions");
  nbc_params->validateParameters(*vp, 0);

  /* create a dummy layout */
  RCP<Layouts> dl = rcp(new Layouts());

  /* temporary variable */
  RCP<PHX::Evaluator<GoalTraits> > ev;

  /* traction bcs */
  if (nbc_params->isSublist("tractions")) {
    RCP<ParameterList> p = rcp(new ParameterList);
    RCP<const ParameterList> tp = rcpFromRef(nbc_params->sublist("tractions"));
    p->set<RCP<Layouts> >("Layouts", dl);
    p->set<RCP<const ParameterList> >("Parameters", tp);
    p->set<RCP<Mesh> >("Mesh", mesh);
    p->set<RCP<Mechanics> >("Mechanics", rcpFromRef(*this));
    ev = rcp(new BCTraction<EvalT, GoalTraits>(*p));
    fm->template registerEvaluator<EvalT>(ev);
    PHX::Tag<typename EvalT::ScalarT> tag("Traction BCs", dl->dummy);
    fm->requireField<EvalT>(tag);
  }
}

}

/* ETI */

template void goal::Mechanics::
register_neumann<goal::GoalTraits::Forward>(FieldManager fm);

template void goal::Mechanics::
register_neumann<goal::GoalTraits::Derivative>(FieldManager fm);
