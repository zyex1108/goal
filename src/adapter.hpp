#ifndef goal_adapter_hpp
#define goal_adapter_hpp

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

class Adapter
{

  public:

    Adapter(
        RCP<const ParameterList> p,
        RCP<Mesh> m,
        RCP<Mechanics> mech,
        RCP<SolutionInfo> s);

    void adapt(unsigned step_number);

  private:

    RCP<const ParameterList> params;
    RCP<Mesh> mesh;
    RCP<Mechanics> mechanics;
    RCP<SolutionInfo> sol_info;

};

RCP<Adapter> adapter_create(
    RCP<const ParameterList> p,
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s);

}

#endif
