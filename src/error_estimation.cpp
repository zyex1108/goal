#include "error_estimation.hpp"
#include "linear_solver.hpp"
#include "mesh.hpp"
#include "mechanics.hpp"
#include "solution_info.hpp"
#include "workset.hpp"
#include "assert_param.hpp"
#include "control.hpp"

namespace goal {

ErrorEstimation::ErrorEstimation(
    RCP<const ParameterList> p,
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s) :
  params(p),
  mesh(m),
  mechanics(mech),
  sol_info(s)
{
}

void ErrorEstimation::set_time(double current, double previous)
{
  t_new = current;
  t_old = previous;
}

void ErrorEstimation::localize()
{
  double t0 = time();
  double t1 = time();
  print("error localized in %f seconds", t1-t0);
}

RCP<ErrorEstimation> error_create(
    RCP<const ParameterList> p,
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s)
{
  RCP<const ParameterList> ep = rcpFromRef(p->sublist("error estimation"));
  return rcp(new ErrorEstimation(ep, m, mech, s));
}

}
