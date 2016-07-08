#include "primal_problem.hpp"
#include "mesh.hpp"
#include "mechanics.hpp"
#include "solution_info.hpp"
#include "workset.hpp"
#include "assert_param.hpp"
#include "control.hpp"

namespace goal {

PrimalProblem::PrimalProblem(
    RCP<const ParameterList> p,
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s) :
  params(p),
  mesh(m),
  mechanics(mech),
  sol_info(s),
  t_new(0.0),
  t_old(0.0)
{
}

void PrimalProblem::set_time(double current, double previous)
{
  t_new = current;
  t_old = previous;
}

void PrimalProblem::compute_residual()
{
}

void PrimalProblem::compute_jacobian()
{
}

void PrimalProblem::solve()
{
}

RCP<PrimalProblem> primal_create(
    RCP<const ParameterList> p,
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s)
{
  RCP<const ParameterList> sp = rcpFromRef(p->sublist("linear algebra"));
  return rcp(new PrimalProblem(sp, m, mech, s));
}

}
