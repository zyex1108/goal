#include "solver_continuation.hpp"
#include "mesh.hpp"
#include "mechanics.hpp"
#include "solution_info.hpp"
#include "primal_problem.hpp"
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
  p->sublist("mesh");
  p->sublist("mechanics");
  p->sublist("linear algebra");
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
  p->validateParameters(*get_valid_params(), 0);
}

SolverContinuation::SolverContinuation(RCP<const ParameterList> p) :
  params(p)
{
  print("--- continuation solver ---");
  validate_params(params);
  bool enable_dynamics = false;
  mesh = mesh_create(params);
  mechanics = mechanics_create(params, mesh, enable_dynamics);
  sol_info = sol_info_create(mesh, enable_dynamics);
  primal = primal_create(params, mesh, mechanics, sol_info);
}

void SolverContinuation::solve()
{
  primal->solve();
}

}
