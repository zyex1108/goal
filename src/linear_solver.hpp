#ifndef goal_linear_solver_hpp
#define goal_linear_solver_hpp

#include "data_types.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

void solve_linear_system(
    RCP<const ParameterList> p,
    RCP<Matrix> A,
    RCP<Vector> x,
    RCP<Vector> b);

}

#endif
