#ifndef goal_solver_goal_continuation_hpp
#define goal_solver_goal_continuation_hpp

#include "solver.hpp"

namespace goal {

class Mesh;
class Mechanics;
class SolutionInfo;
class PrimalProblem;
class DualProblem;
class Output;

class SolverGoalContinuation : public Solver
{
  public:
    SolverGoalContinuation(RCP<const ParameterList> p);
    void solve();
  private:
    RCP<const ParameterList> params;
    RCP<Mesh> mesh;
    RCP<Mechanics> mechanics;
    RCP<SolutionInfo> sol_info;
    RCP<PrimalProblem> primal;
    RCP<DualProblem> dual;
    RCP<Output> output;
    double t_old;
    double t_new;
    double dt;
    unsigned num_steps;
};

}

#endif
