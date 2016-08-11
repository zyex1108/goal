#include "solver_goal_continuation.hpp"
#include "mesh.hpp"
#include "mechanics.hpp"
#include "solution_info.hpp"
#include "primal_problem.hpp"
#include "dual_problem.hpp"
#include "output.hpp"
#include "control.hpp"
#include "assert_param.hpp"

namespace goal {

static RCP<ParameterList> get_valid_params()
{
  RCP<ParameterList> p = rcp(new ParameterList);
  p->set<std::string>("solver type", "");
  p->set<double>("initial time", 0.0);
  p->set<double>("step size", 0.0);
  p->set<unsigned>("num steps", 0.0);
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

SolverGoalContinuation::SolverGoalContinuation(
    RCP<const ParameterList> p) :
  params(p),
  t_old(0.0),
  t_new(0.0),
  dt(0.0),
  num_steps(0)
{
  print("--- goal-oriented adaptive continuation solver ---");
  validate_params(params);
  bool enable_dynamics = false;
  mesh = mesh_create(params);
  mechanics = mechanics_create(params, mesh, enable_dynamics);
  sol_info = sol_info_create(mesh, enable_dynamics);
  primal = primal_create(params, mesh, mechanics, sol_info);
  dual = dual_create(params, mesh, mechanics, sol_info);
  output = output_create(params, mesh, mechanics, sol_info);
  primal->set_coeffs(0.0, 0.0, 1.0);
  dual->set_coeffs(0.0, 0.0, 1.0);
  t_old = params->get<double>("initial time");
  dt = params->get<double>("step size");
  num_steps = params->get<unsigned>("num steps");
  t_new = t_old + dt;
}

void SolverGoalContinuation::solve()
{
  sol_info->ovlp_solution->putScalar(0.0);

  mechanics->build_primal();
  primal->set_time(t_new, t_old);
  primal->solve();

  mechanics->build_dual();
  sol_info->create_dual_vectors(mesh);
  dual->set_time(t_new, t_old);
  dual->solve();

  output->write(0.0);

  sol_info->destroy_dual_vectors();
}

}
