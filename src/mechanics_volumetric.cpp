#include "mechanics.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "mesh.hpp"
#include "control.hpp"
#include "assert_param.hpp"

#include "ev_basis_functions.hpp"
#include "ev_gather_solution.hpp"
#include "ev_dof_interpolation.hpp"
#include "ev_kinematics.hpp"
#include "ev_first_pk.hpp"
#include "ev_elem_size.hpp"
#include "ev_mechanics_residual.hpp"
#include "ev_pressure_residual.hpp"
#include "ev_scatter_residual.hpp"

namespace goal {

template <typename EvalT>
void register_solutions(
    RCP<Layouts> dl,
    RCP<Mesh> mesh,
    Teuchos::Array<std::string>& names,
    unsigned index,
    RCP<PHX::Evaluator<GoalTraits> > ev,
    FieldManager fm)
{
  { /* gather the vector */
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<RCP<Layouts> >("Layouts", dl);
    p->set<RCP<Mesh> >("Mesh", mesh);
    p->set<Teuchos::Array<std::string> >("Sol Names", names);
    p->set<unsigned>("Sol Index", index);
    ev = rcp(new GatherSolution<EvalT, GoalTraits>(*p));
    fm->template registerEvaluator<EvalT>(ev);
  }

  { /* interpolate the degrees of freedom to qps */
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<RCP<Layouts> >("Layouts", dl);
    p->set<Teuchos::Array<std::string> >("DOF Names", names);
    p->set<unsigned>("Sol Index", index);
    p->set<std::string>("BF Name", "BF");
    ev = rcp(new DOFInterpolation<EvalT, GoalTraits>(*p));
    fm->template registerEvaluator<EvalT>(ev);
  }
}

}

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

  /* gather and interpolate the solution vector */
  register_solutions<EvalT>(dl, mesh, var_names[0], 0, ev, fm);

  if (enable_dynamics) {
    /* gather and interpolate the time derivative of the solution */
    register_solutions<EvalT>(dl, mesh, var_names[1], 1, ev, fm);

    /* gather and interpolate the 2nd time derivative of the solution */
    register_solutions<EvalT>(dl, mesh, var_names[2], 2, ev, fm);
  }

  { /* kinematic quantities */
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<RCP<Layouts> >("Layouts", dl);
    p->set<Teuchos::Array<std::string> >("Disp Names", fields["disp"]);
    p->set<std::string>("Def Grad Name", "F");
    p->set<std::string>("Det Def Grad Name", "J");
    ev = rcp(new Kinematics<EvalT, GoalTraits>(*p));
    fm->template registerEvaluator<EvalT>(ev);
  }

  /* material parameters */
  RCP<const ParameterList> mp = rcpFromRef(params->sublist(set));

  /* temperature parameters */
  RCP<const ParameterList> tp;
  if (have_temperature) tp = rcpFromRef(params->sublist("temperature"));

  this->template register_model<EvalT>(set, mp, tp, fm);

  { /* first piola kirchhoff stress tensor */
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<RCP<Layouts> >("Layouts", dl);
    p->set<bool>("Small Strain", small_strain);
    p->set<bool>("Have Pressure", have_pressure_eq);
    p->set<std::string>("Def Grad Name", "F");
    p->set<std::string>("Det Def Grad Name", "J");
    p->set<std::string>("Cauchy Name", "cauchy");
    if (have_pressure_eq)
      p->set<std::string>("Pressure Name", "p");
    p->set<std::string>("First PK Name", "first_pk");
    ev = rcp(new FirstPK<EvalT, GoalTraits>(*p));
    fm->template registerEvaluator<EvalT>(ev);
  }

  { /* mechanics residual */
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<RCP<Layouts> >("Layouts", dl);
    p->set<RCP<Mesh> >("Mesh", mesh);
    p->set<std::string>("BF Name", "BF");
    p->set<Teuchos::Array<std::string> >("Disp Names", fields["disp"]);
    p->set<std::string>("Weighted Dv Name", "wDv");
    p->set<std::string>("First PK Name", "first_pk");
    p->set<bool>("Enable Dynamics", enable_dynamics);
    if (enable_dynamics) {
      p->set<double>("Density", mp->get<double>("rho"));
      p->set<Teuchos::Array<std::string> >("Acc Names", fields["acc"]);
    }
    p->set<bool>("Have Body Force", have_body_force);
    if (have_body_force) {
      Teuchos::Array<std::string> body_force;
      body_force = params->get<Teuchos::Array<std::string> >("body force");
      p->set<Teuchos::Array<std::string> >("Body Force", body_force);
      p->set<double>("Density", mp->get<double>("rho"));
    }
    ev = rcp(new MechanicsResidual<EvalT, GoalTraits>(*p));
    fm->template registerEvaluator<EvalT>(ev);
  }

  if (have_pressure_eq) { /* element size */
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<RCP<Layouts> >("Layouts", dl);
    p->set<std::string>("BF Name", "BF");
    p->set<std::string>("Size Name", "size");
    ev = rcp(new ElemSize<EvalT, GoalTraits>(*p));
    fm->template registerEvaluator<EvalT>(ev);
  }

  if (have_pressure_eq) { /* pressure residual */
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<RCP<Layouts> >("Layouts", dl);
    p->set<RCP<const ParameterList> >("Material Params", mp);
    p->set<bool>("Small Strain", small_strain);
    p->set<std::string>("Pressure Name", "p");
    p->set<std::string>("Size Name", "size");
    p->set<std::string>("Weighted Dv Name", "wDv");
    p->set<std::string>("BF Name", "BF");
    p->set<std::string>("Def Grad Name", "F");
    p->set<std::string>("Det Def Grad Name", "J");
    p->set<std::string>("Cauchy Name", "cauchy");
    ev = rcp(new PressureResidual<EvalT, GoalTraits>(*p));
    fm->template registerEvaluator<EvalT>(ev);
  }

  { /* scatter residuals */
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<RCP<Layouts> >("Layouts", dl);
    p->set<RCP<Mesh> >("Mesh", mesh);
    p->set<Teuchos::Array<std::string> >("DOF Names", var_names[0]);
    ev = rcp(new ScatterResidual<EvalT, GoalTraits>(*p));
    fm->template registerEvaluator<EvalT>(ev);
    PHX::Tag<typename EvalT::ScalarT> tag("Scatter Residual", dl->dummy);
    fm->requireField<EvalT>(tag);
  }

  if (is_dual)
    this->template register_qoi<EvalT>(set, fm);

}

/* ETI */

template void goal::Mechanics::
register_volumetric<goal::GoalTraits::Forward>(
    std::string const& set, FieldManager fm);

template void goal::Mechanics::
register_volumetric<goal::GoalTraits::Derivative>(
    std::string const& set, FieldManager fm);
