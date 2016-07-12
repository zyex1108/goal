#ifndef goal_solver_hpp
#define goal_solver_hpp

#include <Teuchos_RCP.hpp>

namespace Teuchos {
class ParameterList;
}

namespace goal {

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::ParameterList;

class Solver
{
  public:
    virtual ~Solver();
    virtual void solve() = 0;
};

RCP<Solver> solver_create(RCP<const ParameterList> p);

}

#endif
