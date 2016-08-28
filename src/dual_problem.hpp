#ifndef goal_dual_problem_hpp
#define goal_dual_problem_hpp

#include <Teuchos_RCP.hpp>

namespace Teuchos {
class ParameterList;
}

namespace goal {

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::ParameterList;

class Mesh;
class Mechanics;
class SolutionInfo;

class DualProblem
{
  public:

    DualProblem(
        RCP<const ParameterList> p,
        RCP<Mesh> m,
        RCP<Mechanics> mech,
        RCP<SolutionInfo> s);

    void set_time(double current, double previous);

    void set_coeffs(double alpha, double beta, double gamma);

    void solve();

  private:

    RCP<const ParameterList> params;
    RCP<Mesh> mesh;
    RCP<Mechanics> mechanics;
    RCP<SolutionInfo> sol_info;

    double t_new;
    double t_old;

    double alpha;
    double beta;
    double gamma;

    void compute_jacobian();

};

RCP<DualProblem> dual_create(
    RCP<const ParameterList> p,
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s);

}

#endif
