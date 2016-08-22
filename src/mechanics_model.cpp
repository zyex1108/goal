#include "mechanics.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "mesh.hpp"
#include "control.hpp"
#include "assert_param.hpp"

#include "ev_model_elastic.hpp"
#include "ev_model_j2.hpp"
#include "ev_model_creep.hpp"

template <typename EvalT>
void goal::Mechanics::register_model(
    std::string const& set,
    RCP<const ParameterList> material_params,
    RCP<const ParameterList> temperature_params,
    FieldManager fm)
{
  /* do some work to create a data layout */
  unsigned ws_size = mesh->get_ws_size();
  unsigned n_nodes = mesh->get_num_elem_nodes();
  unsigned n_qps = mesh->get_num_elem_qps();
  unsigned n_dims = mesh->get_num_dims();
  RCP<Layouts> dl = rcp(new Layouts(ws_size, n_nodes, n_qps, n_dims));

  /* temporary variable */
  RCP<PHX::Evaluator<GoalTraits> > ev;

  if (model == "linear elastic") { /* elastic model */
    small_strain = true;
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<RCP<Layouts> >("Layouts", dl);
    p->set<RCP<StateFields> >("State Fields", state_fields);
    p->set<RCP<const ParameterList> >("Material Params",material_params);
    p->set<RCP<const ParameterList> >("Temperature Params",temperature_params);
    p->set<Teuchos::Array<std::string> >("Disp Names", fields["disp"]);
    p->set<std::string>("Cauchy Name", "cauchy");
    ev = rcp(new ModelElastic<EvalT, GoalTraits>(*p));
    fm->template registerEvaluator<EvalT>(ev);
  }

  else if (model == "j2") { /* j2 model */
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<RCP<Layouts> >("Layouts", dl);
    p->set<RCP<StateFields> >("State Fields", state_fields);
    p->set<RCP<const ParameterList> >("Material Params",material_params);
    p->set<RCP<const ParameterList> >("Temperature Params",temperature_params);
    p->set<std::string>("Def Grad Name", "F");
    p->set<std::string>("Det Def Grad Name", "J");
    p->set<std::string>("Cauchy Name", "cauchy");
    ev = rcp(new ModelJ2<EvalT, GoalTraits>(*p));
    fm->template registerEvaluator<EvalT>(ev);
  }

  else if (model == "creep") { /* creep model */
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<RCP<Layouts> >("Layouts", dl);
    p->set<RCP<StateFields> >("State Fields", state_fields);
    p->set<std::string>("Def Grad Name", "F");
    p->set<std::string>("Det Def Grad Name", "J");
    p->set<std::string>("Cauchy Name", "cauchy");
    p->set<RCP<const ParameterList> >("Material Params",material_params);
    ev = rcp(new ModelCreep<EvalT, GoalTraits>(*p));
    fm->template registerEvaluator<EvalT>(ev);
  }

  else
    fail("unknown material model: %s", model.c_str());
}

/* ETI */

template void goal::Mechanics::
register_model<goal::GoalTraits::Forward>(
    std::string const& set,
    RCP<const ParameterList> material_params,
    RCP<const ParameterList> temperature_params,
    FieldManager fm);

template void goal::Mechanics::
register_model<goal::GoalTraits::Derivative>(
    std::string const& set,
    RCP<const ParameterList> material_params,
    RCP<const ParameterList> temperature_params,
    FieldManager fm);
