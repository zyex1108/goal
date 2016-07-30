#include "primal_problem.hpp"
#include "linear_solver.hpp"
#include "mesh.hpp"
#include "mechanics.hpp"
#include "solution_info.hpp"
#include "workset.hpp"
#include "assert_param.hpp"
#include "control.hpp"

namespace goal {

static RCP<ParameterList> get_valid_params()
{
  RCP<ParameterList> p = rcp(new ParameterList);
  p->set<double>("linear: tolerance", 0.0);
  p->set<unsigned>("linear: max iters", 0);
  p->set<unsigned>("linear: krylov size", 0);
  p->set<double>("nonlinear: tolerance", 0.0);
  p->set<unsigned>("nonlinear: max iters", 0);
  return p;
}

static void validate_params(RCP<const ParameterList> p)
{
  assert_param(p, "linear: tolerance");
  assert_param(p, "linear: max iters");
  assert_param(p, "linear: krylov size");
  assert_param(p, "nonlinear: tolerance");
  assert_param(p, "nonlinear: max iters");
  p->validateParameters(*get_valid_params(), 0);
}

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
  t_old(0.0),
  alpha(0.0),
  beta(0.0),
  gamma(0.0)
{
  validate_params(params);
  tolerance = params->get<double>("nonlinear: tolerance");
  max_iters = params->get<unsigned>("nonlinear: max iters");
}

struct PrimalInfo
{
  double t_new;
  double t_old;
  double alpha;
  double beta;
  double gamma;
};

void PrimalProblem::set_time(double current, double previous)
{
  t_new = current;
  t_old = previous;
}

void PrimalProblem::set_coeffs(double a, double b, double c)
{
  alpha = a;
  beta = b;
  gamma = c;
}

static void load_overlap_solution(Workset& ws, RCP<SolutionInfo> s)
{
  ws.u = s->ovlp_solution;
  ws.r = s->ovlp_residual;
  ws.J = s->ovlp_jacobian;
}

static void load_owned_solution(Workset& ws, RCP<SolutionInfo> s)
{
  ws.u = s->owned_solution;
  ws.r = s->owned_residual;
  ws.J = s->owned_jacobian;
}

static void load_primal_info(Workset& ws, PrimalInfo* info)
{
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

static void compute_volumetric_residual(
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s,
    PrimalInfo* info)
{
  typedef GoalTraits::Residual R;
  FieldManagers f = mech->get_volumetric();
  Workset ws;
  load_overlap_solution(ws, s);
  load_primal_info(ws, info);
  for (unsigned set_idx=0; set_idx < f.size(); ++set_idx) {
    f[set_idx]->postRegistrationSetupForType<R>(NULL);
    unsigned num_ws = m->get_num_worksets(set_idx);
    for (unsigned ws_idx=0; ws_idx < num_ws; ++ws_idx) {
      load_mesh_info(ws, m, set_idx, ws_idx);
      f[set_idx]->evaluateFields<R>(ws);
    }
  }
}

static void compute_neumann_residual(
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s,
    PrimalInfo* info)
{
  typedef GoalTraits::Residual R;
  FieldManager f = mech->get_neumann();
  Workset ws;
  load_overlap_solution(ws, s);
  load_primal_info(ws, info);
  f->postRegistrationSetupForType<R>(NULL);
  f->evaluateFields<R>(ws);
}

static void compute_dirichlet_residual(
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s,
    PrimalInfo* info)
{
  typedef GoalTraits::Residual R;
  FieldManager f = mech->get_dirichlet();
  Workset ws;
  load_owned_solution(ws, s);
  load_primal_info(ws, info);
  f->postRegistrationSetupForType<R>(NULL);
  f->evaluateFields<R>(ws);
}

void PrimalProblem::compute_residual()
{
  double t0 = time();
  sol_info->scatter_solution();
  sol_info->owned_residual->putScalar(0.0);
  sol_info->ovlp_residual->putScalar(0.0);
  PrimalInfo primal_info = {t_new,t_old,alpha,beta,gamma};
  compute_volumetric_residual(mesh, mechanics, sol_info, &primal_info);
  compute_neumann_residual(mesh, mechanics, sol_info, &primal_info);
  sol_info->gather_residual();
  compute_dirichlet_residual(mesh, mechanics, sol_info, &primal_info);
  double t1 = time();
  print("  residual computed in %f seconds", t1-t0);
}

static void compute_volumetric_jacobian(
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s,
    PrimalInfo* info)
{
  typedef GoalTraits::Jacobian J;
  FieldManagers f = mech->get_volumetric();
  Workset ws;
  load_overlap_solution(ws, s);
  load_primal_info(ws, info);
  for (unsigned set_idx=0; set_idx < f.size(); ++set_idx) {
    std::vector<PHX::index_size_type> dd;
    dd.push_back(m->get_num_elem_dofs());
    f[set_idx]->setKokkosExtendedDataTypeDimensions<J>(dd);
    f[set_idx]->postRegistrationSetupForType<J>(NULL);
    unsigned num_ws = m->get_num_worksets(set_idx);
    for (unsigned ws_idx=0; ws_idx < num_ws; ++ws_idx) {
      load_mesh_info(ws, m, set_idx, ws_idx);
      f[set_idx]->evaluateFields<J>(ws);
    }
  }
}

static void compute_neumann_jacobian(
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s,
    PrimalInfo* info)
{
  typedef GoalTraits::Jacobian J;
  FieldManager f = mech->get_neumann();
  Workset ws;
  load_overlap_solution(ws, s);
  load_primal_info(ws, info);
  std::vector<PHX::index_size_type> dd;
  dd.push_back(m->get_num_elem_dofs());
  f->setKokkosExtendedDataTypeDimensions<J>(dd);
  f->postRegistrationSetupForType<J>(NULL);
  f->evaluateFields<J>(ws);
}

static void compute_dirichlet_jacobian(
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s,
    PrimalInfo* info)
{
  typedef GoalTraits::Jacobian J;
  FieldManager f = mech->get_dirichlet();
  Workset ws;
  load_owned_solution(ws, s);
  load_primal_info(ws, info);
  std::vector<PHX::index_size_type> dd;
  dd.push_back(m->get_num_elem_dofs());
  f->setKokkosExtendedDataTypeDimensions<J>(dd);
  f->postRegistrationSetupForType<J>(NULL);
  f->evaluateFields<J>(ws);
}

void PrimalProblem::compute_jacobian()
{
  double t0 = time();
  sol_info->scatter_solution();
  sol_info->owned_residual->putScalar(0.0);
  sol_info->ovlp_residual->putScalar(0.0);
  sol_info->owned_jacobian->resumeFill();
  sol_info->ovlp_jacobian->resumeFill();
  sol_info->owned_jacobian->setAllToScalar(0.0);
  sol_info->ovlp_jacobian->setAllToScalar(0.0);
  PrimalInfo primal_info = {t_new,t_old,alpha,beta,gamma};
  compute_volumetric_jacobian(mesh, mechanics, sol_info, &primal_info);
  compute_neumann_jacobian(mesh, mechanics, sol_info, &primal_info);
  sol_info->ovlp_jacobian->fillComplete();
  sol_info->gather_residual();
  sol_info->gather_jacobian();
  compute_dirichlet_jacobian(mesh, mechanics, sol_info, &primal_info);
  sol_info->owned_jacobian->fillComplete();
  double t1 = time();
  print("  jacobian computed in %f seconds", t1-t0);
}

void PrimalProblem::solve()
{
  print("solving primal model");
  RCP<Matrix> J = sol_info->owned_jacobian;
  RCP<Vector> u = sol_info->owned_solution->getVectorNonConst(0);
  RCP<Vector> r = sol_info->owned_residual;
  RCP<Vector> du = rcp(new Vector(mesh->get_owned_map()));
  unsigned iter=1;
  bool converged = false;
  while ((iter <= max_iters) && (! converged)) {
    print(" (%d) newton iteration", iter);
    compute_jacobian();
    r->scale(-1.0);
    du->putScalar(0.0);
    solve_linear_system(params, J, du, r);
    u->update(1.0, *du, 1.0);
    compute_residual();
    double norm = r->norm2();
    print("  ||r|| = %e", norm);
    if (norm < tolerance) converged = true;
    iter++;
  }
  if ((iter > max_iters) && (!converged))
    fail("newton's method failed in %u iterations", max_iters);
  du = Teuchos::null;
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
