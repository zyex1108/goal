#include "solver_trapezoid.hpp"
#include "mesh.hpp"
#include "mechanics.hpp"
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

SolverTrapezoid::SolverTrapezoid(RCP<const ParameterList> p) :
  params(p)
{
  print("--- trapezoid solver ---");
  validate_params(params);
  mesh = mesh_create(params);
  mechanics = mechanics_create(params, mesh, true);
}

void SolverTrapezoid::solve()
{
}

}
