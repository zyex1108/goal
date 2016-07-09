#ifndef goal_solver_trapezoid_hpp
#define goal_solver_trapezoid_hpp

#include "solver.hpp"

namespace goal {

class Mesh;
class Mechanics;
class SolutionInfo;
class PrimalProblem;

class SolverTrapezoid : public Solver
{
  public:
    SolverTrapezoid(RCP<const ParameterList> p);
    void solve();
  private:
    RCP<const ParameterList> params;
    RCP<Mesh> mesh;
    RCP<Mechanics> mechanics;
    RCP<SolutionInfo> sol_info;
    RCP<PrimalProblem> primal;
};

}

#endif
