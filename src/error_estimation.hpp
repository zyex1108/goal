#ifndef goal_error_estimation_hpp
#define goal_error_estimation_hpp

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

class ErrorEstimation
{
  public:

    ErrorEstimation(
        RCP<const ParameterList> p,
        RCP<Mesh> m,
        RCP<Mechanics> mech,
        RCP<SolutionInfo> s);

    void set_time(double current, double previous);

    void localize();

  private:

    RCP<const ParameterList> params;
    RCP<Mesh> mesh;
    RCP<Mechanics> mechanics;
    RCP<SolutionInfo> sol_info;

    double t_new;
    double t_old;

    void compute_error();

};

RCP<ErrorEstimation> error_create(
    RCP<const ParameterList> p,
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s);

}

#endif
