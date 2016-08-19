#include "solver_trapezoid.hpp"
#include "mesh.hpp"
#include "mechanics.hpp"
#include "solution_info.hpp"
#include "initial_condition.hpp"
#include "primal_problem.hpp"
#include "output.hpp"
#include "assert_param.hpp"
#include "control.hpp"

#include <Teuchos_ParameterList.hpp>

namespace goal {

static RCP<ParameterList> get_valid_params()
{
  RCP<ParameterList> p = rcp(new ParameterList);
  Teuchos::Array<std::string> dummy(0);
  p->set<std::string>("solver type", "");
  p->set<double>("initial time", 0.0);
  p->set<double>("step size", 0.0);
  p->set<unsigned>("num steps", 0);
  p->set<Teuchos::Array<std::string> >("ic", dummy);
  p->set<Teuchos::Array<std::string> >("ic dot", dummy);
  p->sublist("mesh");
  p->sublist("mechanics");
  p->sublist("linear algebra");
  p->sublist("output");
  return p;
}

static void validate_params(RCP<const ParameterList> p)
{
  assert_param(p, "solver type");
  assert_param(p, "initial time");
  assert_param(p, "step size");
  assert_param(p, "num steps");
  assert_sublist(p, "mesh");
  assert_sublist(p, "mechanics");
  assert_sublist(p, "output");
  p->validateParameters(*get_valid_params(), 0);
}

SolverTrapezoid::SolverTrapezoid(RCP<const ParameterList> p) :
  params(p),
  t_old(0.0),
  t_new(0.0),
  dt(0.0),
  num_steps(0)
{
  print("--- trapezoid solver ---");
  validate_params(params);
  bool enable_dynamics = true;
  mesh = mesh_create(params);
  mechanics = mechanics_create(params, mesh, enable_dynamics);
  sol_info = sol_info_create(mesh, enable_dynamics);
  primal = primal_create(params, mesh, mechanics, sol_info);
  output = output_create(params, mesh, mechanics, sol_info);
  set_initial_conditions(params, mesh, mechanics, sol_info);
  t_old = params->get<double>("initial time");
  dt = params->get<double>("step size");
  num_steps = params->get<unsigned>("num steps");
  t_new = t_old + dt;
}

void SolverTrapezoid::solve()
{
  mechanics->build_primal();

  /* create the predictor vectors */
  RCP<const Map> map = mesh->get_owned_map();
  x_pred = rcp(new Vector(map));
  v_pred = rcp(new Vector(map));

  /* time loop */
  for (unsigned step=1; step <= num_steps; ++step) {

    print("*** Time Step: (%u)", step);
    print("*** from time: %f", t_old);
    print("*** to time:   %f", t_new);

    /* compute fad coefficients */
    double alpha = 4.0/(dt*dt);
    double beta = 2.0/dt;
    double gamma = 1.0;

    /* get the solution vectors */
    RCP<Vector> x = sol_info->owned_solution->getVectorNonConst(0);
    RCP<Vector> v = sol_info->owned_solution->getVectorNonConst(1);
    RCP<Vector> a = sol_info->owned_solution->getVectorNonConst(2);

    /* predictor phase */
    x_pred->assign(*x);
    x_pred->update(dt, *v, alpha, *a, 1.0);
    v_pred->assign(*v);
    v_pred->update(beta, *v, 1.0);
    a->assign(*v_pred);

    /* solve the primal model */
    primal->set_coeffs(alpha, beta, gamma);
    primal->set_time(t_new, t_old);
    primal->solve();

    /* corrector phase */
    a->update(alpha, *x, -alpha, *x_pred, 0.0);
    v->update(beta, *x, -beta, *v_pred, 0.0);

    /* write output */
    output->write(t_new);

    /* updates */
    t_old = t_new;
    t_new = t_new + dt;
    mechanics->update_state();

  }

}

}
