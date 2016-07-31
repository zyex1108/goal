#include "solver.hpp"
#include "solver_continuation.hpp"
#include "solver_goal_continuation.hpp"
#include "assert_param.hpp"
#include "control.hpp"

#include <Teuchos_ParameterList.hpp>

namespace goal {

Solver::~Solver()
{
}

RCP<Solver> solver_create(RCP<const ParameterList> p)
{
  assert_param(p, "solver type");
  std::string const& type = p->get<std::string>("solver type");
  Teuchos::RCP<Solver> solver;
  if (type == "continuation")
    solver = rcp(new SolverContinuation(p));
  else if (type == "goal-oriented continuation")
    solver = rcp(new SolverGoalContinuation(p));
  else
    fail("unknown solver type: %s", type.c_str());
  return solver;
}

}
