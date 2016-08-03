#include "linear_solver.hpp"
#include "control.hpp"

#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Ifpack2_Factory.hpp>

namespace goal {

typedef Tpetra::MultiVector<ST, LO, GO, KNode> MV;
typedef Tpetra::Operator<ST, LO, GO, KNode> OP;
typedef Tpetra::RowMatrix<ST, LO, GO, KNode> RM;
typedef Belos::LinearProblem<ST, MV, OP> LinearProblem;
typedef Belos::SolverManager<ST, MV, OP> Solver;
typedef Belos::BlockGmresSolMgr<ST, MV, OP> GmresSolver;
typedef Tpetra::Operator<ST, LO, GO, KNode> Prec;
typedef Ifpack2::Preconditioner<ST, LO, GO, KNode> IfpackPrec;

static RCP<ParameterList> get_ifpack2_params()
{
  RCP<ParameterList> p = rcp(new ParameterList);
  p->set("fact: drop tolerance", 0.0);
  p->set("fact: ilut level-of-fill", 1.0);
  return p;
}

static RCP<ParameterList> get_belos_params(RCP<const ParameterList> in)
{
  RCP<ParameterList> p = rcp(new ParameterList);
  int max_iters = in->get<unsigned>("linear: max iters");
  int krylov = in->get<unsigned>("linear: krylov size");
  double tol = in->get<double>("linear: tolerance");
  p->set("Block Size", 1);
  p->set("Num Blocks", krylov);
  p->set("Maximum Iterations", max_iters);
  p->set("Convergence Tolerance", tol);
  p->set("Orthogonalization", "DGKS");
  return p;
}

static RCP<Prec> build_ifpack2_prec(RCP<Matrix> A)
{
  RCP<ParameterList> p = get_ifpack2_params();
  Ifpack2::Factory factory;
  RCP<IfpackPrec> prec = factory.create<RM>("ILUT", A);
  prec->setParameters(*p);
  prec->initialize();
  prec->compute();
  return prec;
}

static RCP<Prec> build_precond(RCP<const ParameterList> p, RCP<Matrix> A)
{
  (void)(p);
  return build_ifpack2_prec(A);
}

static RCP<Solver> build_solver(
    RCP<const ParameterList> in,
    RCP<Prec> P,
    RCP<Matrix> A,
    RCP<Vector> x,
    RCP<Vector> b)
{
  RCP<ParameterList> p = get_belos_params(in);
  RCP<LinearProblem> problem = rcp(new LinearProblem(A,x,b));
  problem->setLeftPrec(P);
  problem->setProblem();
  RCP<Solver> solver = rcp(new GmresSolver(problem, p));
  return solver;
}

void solve_linear_system(
    RCP<const ParameterList> in,
    RCP<Matrix> A,
    RCP<Vector> x,
    RCP<Vector> b)
{
  double t0 = time();
  RCP<Prec> P = build_precond(in, A);
  RCP<Solver> solver = build_solver(in, P, A, x, b);
  solver->solve();
  unsigned iters = solver->getNumIters();
  double t1 = time();
  if (iters >= in->get<unsigned>("linear: max iters"))
    print("  linear solve failed to converge in %d iterations\n"
          "  continuing using the incomplete solve...", iters);
  else
    print("  linear system solved in %d iterations", iters);
  print("  linear system solved in %f seconds", t1-t0);
}

}
