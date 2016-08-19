#ifndef goal_solver_trapezoid_hpp
#define goal_solver_trapezoid_hpp

#include "solver.hpp"
#include "data_types.hpp"

namespace goal {

class Mesh;
class Mechanics;
class SolutionInfo;
class PrimalProblem;
class Output;

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
    RCP<Output> output;
    RCP<Vector> x_pred;
    RCP<Vector> v_pred;
    double t_old;
    double t_new;
    double dt;
    unsigned num_steps;
};

}

#endif
