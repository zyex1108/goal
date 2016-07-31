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

struct DualInfo
{
  double t_new;
  double t_old;
  double alpha;
  double beta;
  double gamma;
};

void DualProblem::set_time(double current, double previous)
{
  t_new = current;
  t_old = previous;
}

void DualProblem::set_coeffs(double a, double b, double c)
{
  alpha = a;
  beta = b;
  gamma = c;
}

static void load_overlap_solution(Workset& ws, RCP<SolutionInfo> s)
{
  ws.u = s->ovlp_solution;
  ws.q = s->ovlp_qoi;
  ws.J = s->ovlp_jacobian;
}

static void load_owned_solution(Workset& ws, RCP<SolutionInfo> s)
{
  ws.u = s->owned_solution;
  ws.q = s->owned_qoi;
  ws.J = s->owned_jacobian;
}

static void load_dual_info(Workset& ws, DualInfo* info)
{
  ws.is_adjoint = true;
  ws.t_old = info->t_old;
  ws.t_new = info->t_new;
  ws.alpha = info->alpha;
  ws.beta = info->beta;
  ws.gamma = info->gamma;
}

static void load_mesh_info(
    Workset& ws,
    RCP<Mesh> m,
    const unsigned set_idx,
    const unsigned ws_idx)
{
  std::string const& set = m->get_elem_set_name(set_idx);
  ws.set = set;
  ws.ents = m->get_elems(set, ws_idx);
  ws.size = ws.ents.size();
}

static void compute_volumetric_jacobian(
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s,
    DualInfo* info)
{
  typedef GoalTraits::Derivative D;
  FieldManagers f = mech->get_volumetric();
  Workset ws;
  load_overlap_solution(ws, s);
  load_dual_info(ws, info);
  for (unsigned set_idx=0; set_idx < f.size(); ++set_idx) {
    std::vector<PHX::index_size_type> dd;
    dd.push_back(m->get_num_elem_dofs());
    f[set_idx]->setKokkosExtendedDataTypeDimensions<D>(dd);
    f[set_idx]->postRegistrationSetupForType<D>(NULL);
    unsigned num_ws = m->get_num_worksets(set_idx);
    for (unsigned ws_idx=0; ws_idx < num_ws; ++ws_idx) {
      load_mesh_info(ws, m, set_idx, ws_idx);
      f[set_idx]->evaluateFields<D>(ws);
    }
  }
}

static void compute_dirichlet_jacobian(
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s,
    DualInfo* info)
{
  typedef GoalTraits::Derivative D;
  FieldManager f = mech->get_dirichlet();
  Workset ws;
  load_owned_solution(ws, s);
  load_dual_info(ws, info);
  std::vector<PHX::index_size_type> dd;
  dd.push_back(m->get_num_elem_dofs());
  f->setKokkosExtendedDataTypeDimensions<D>(dd);
  f->postRegistrationSetupForType<D>(NULL);
  f->evaluateFields<D>(ws);
}

void DualProblem::compute_jacobian()
{
  double t0 = time();
  sol_info->scatter_solution();
  sol_info->owned_qoi->putScalar(0.0);
  sol_info->ovlp_qoi->putScalar(0.0);
  sol_info->owned_jacobian->resumeFill();
  sol_info->ovlp_jacobian->resumeFill();
  sol_info->owned_jacobian->setAllToScalar(0.0);
  sol_info->ovlp_jacobian->setAllToScalar(0.0);
  DualInfo dual_info = {t_new,t_old,alpha,beta,gamma};
  compute_volumetric_jacobian(mesh, mechanics, sol_info, &dual_info);
  sol_info->ovlp_jacobian->fillComplete();
  sol_info->gather_qoi();
  sol_info->gather_jacobian();
  compute_dirichlet_jacobian(mesh, mechanics, sol_info, &dual_info);
  sol_info->owned_jacobian->fillComplete();
  double t1 = time();
  print("  jacobian transpose computed in %f seconds", t1-t0);
}

void DualProblem::solve()
{
  print("solving dual model");
  compute_jacobian();
  RCP<Vector> z = sol_info->owned_dual;
  RCP<Vector> q = sol_info->owned_qoi;
  RCP<Matrix> J = sol_info->owned_jacobian;
  solve_linear_system(params, J, z, q);
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
