#ifndef goal_initial_condition_hpp
#define goal_initial_condition_hpp

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class Mesh;
class Mechanics;
class SolutionInfo;

void set_initial_conditions(
    RCP<const ParameterList> p,
    RCP<Mesh> mesh,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> info);

}

#endif
