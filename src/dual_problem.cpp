#include "dual_problem.hpp"
#include "linear_solver.hpp"
#include "mesh.hpp"
#include "mechanics.hpp"
#include "solution_info.hpp"
#include "workset.hpp"
#include "assert_param.hpp"
#include "control.hpp"

namespace goal {

static void validate_params(RCP<const ParameterList> p)
{
  assert_param(p, "linear: tolerance");
  assert_param(p, "linear: max iters");
  assert_param(p, "linear: krylov size");
}

DualProblem::DualProblem(
    RCP<const ParameterList> p,
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s) :
  params(p),
  mesh(m),
  mechanics(mech),
  sol_info(s),
  t_new(0.0),
  t_old(0.0),
  alpha(0.0),
  beta(0.0),
  gamma(0.0)
{
  validate_params(params);
}

void DualProblem::solve()
{
}

RCP<DualProblem> dual_create(
    RCP<const ParameterList> p,
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s)
{
  RCP<const ParameterList> dp = rcpFromRef(p->sublist("linear algebra"));
  return rcp(new DualProblem(dp, m, mech, s));
}

}
