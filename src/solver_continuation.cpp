#include "solver_continuation.hpp"
#include "mesh.hpp"
#include "mechanics.hpp"
#include "solution_info.hpp"
#include "primal_problem.hpp"
#include "output.hpp"
#include "assert_param.hpp"
#include "control.hpp"

#include <Teuchos_ParameterList.hpp>

namespace goal {

static RCP<ParameterList> get_valid_params()
{
  RCP<ParameterList> p = rcp(new ParameterList);
  p->set<std::string>("solver type", "");
  p->set<double>("initial time", 0.0);
  p->set<double>("step size", 0.0);
  p->set<unsigned>("num steps", 0);
  p->set<double>("regression: val", 0.0);
  p->set<double>("regression: tol", 0.0);
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
  assert_param(p, "mesh");
  assert_param(p, "mechanics");
  assert_param(p, "output");
  p->validateParameters(*get_valid_params(), 0);
}

SolverContinuation::SolverContinuation(RCP<const ParameterList> p) :
  params(p),
  t_old(0.0),
  t_new(0.0),
  dt(0.0),
  num_steps(0)
{
  print("--- continuation solver ---");
  validate_params(params);
  bool enable_dynamics = false;
  mesh = mesh_create(params);
  mechanics = mechanics_create(params, mesh, enable_dynamics);
  sol_info = sol_info_create(mesh, enable_dynamics);
  primal = primal_create(params, mesh, mechanics, sol_info);
  output = output_create(params, mesh, mechanics, sol_info);
  primal->set_coeffs(0.0, 0.0, 1.0);
  t_old = params->get<double>("initial time");
  dt = params->get<double>("step size");
  num_steps = params->get<unsigned>("num steps");
  t_new = t_old + dt;
}

static void check_regression(
    RCP<const ParameterList> p,
    RCP<SolutionInfo> s)
{
  double tol = p->get<double>("regression: tol");
  double expected = p->get<double>("regression: val");
  RCP<const Vector> x = s->owned_solution->getVector(0);
  double computed = x->meanValue();
  print("expected solution average: %.15f", expected);
  print("computed solution average: %.15f", computed);
  CHECK(std::abs(computed-expected) < tol);
}

void SolverContinuation::solve()
{
  mechanics->build_primal();
  sol_info->ovlp_solution->putScalar(0.0);
  for (unsigned step=1; step <= num_steps; ++step) {
    print("*** Continuation Step: (%u)", step);
    print("*** from time:         %f", t_old);
    print("*** to time:           %f", t_new);
    primal->set_time(t_new, t_old);
    primal->solve();
    output->write(t_new);
    t_old = t_new;
    t_new = t_new + dt;
    mechanics->update_state();
  }
  if (params->isParameter("regression: val"))
    check_regression(params, sol_info);
}

}
