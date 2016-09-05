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
  sol_info(s),
  t_new(0.0),
  t_old(0.0)
{
}

void ErrorEstimation::set_time(double current, double previous)
{
  t_new = current;
  t_old = previous;
}

static void load_time_info(
    Workset& ws,
    const double t_new,
    const double t_old)
{
  ws.t_new = t_new;
  ws.t_old = t_old;
}

static void load_overlap_solution(Workset& ws, RCP<SolutionInfo> s)
{
  ws.u = s->ovlp_solution;
  ws.z = s->ovlp_dual;
  ws.r = s->ovlp_residual;
}

static void load_mesh_info(
    Workset& ws, RCP<Mesh> m,
    const unsigned set_idx,
    const unsigned ws_idx)
{
  std::string const& set = m->get_elem_set_name(set_idx);
  ws.set = set;
  ws.ents = m->get_elems(set, ws_idx);
  ws.size = ws.ents.size();
}

static void compute_volumetric_error(
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s,
    const double t_new,
    const double t_old)
{
  typedef GoalTraits::Forward F;
  FieldManagers f = mech->get_volumetric();
  Workset ws;
  load_overlap_solution(ws, s);
  load_time_info(ws, t_new, t_old);
  for (unsigned set_idx=0; set_idx < f.size(); ++set_idx) {
    f[set_idx]->postRegistrationSetupForType<F>(NULL);
    unsigned num_ws = m->get_num_worksets(set_idx);
    for (unsigned ws_idx=0; ws_idx < num_ws; ++ws_idx) {
      load_mesh_info(ws, m, set_idx, ws_idx);
      f[set_idx]->evaluateFields<F>(ws);
    }
  }
}

void ErrorEstimation::localize()
{
  double t0 = time();
  sol_info->scatter_solution();
  sol_info->owned_residual->putScalar(0.0);
  sol_info->ovlp_residual->putScalar(0.0);
  compute_volumetric_error(mesh, mechanics, sol_info, t_new, t_old);
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
