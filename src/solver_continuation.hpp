#ifndef goal_solver_continuation_hpp
#define goal_solver_continuation_hpp

#include "solver.hpp"

namespace goal {

class Mesh;
class Mechanics;

class SolverContinuation : public Solver
{
  public:
    SolverContinuation(RCP<const ParameterList> p);
    void solve();
  private:
    RCP<const ParameterList> params;
    RCP<Mesh> mesh;
    RCP<Mechanics> mechanics;
};

}

#endif
