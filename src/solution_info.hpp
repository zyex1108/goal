#ifndef goal_solution_info_hpp
#define goal_solution_info_hpp

#include "data_types.hpp"

namespace goal {

class SolutionInfo
{
  public:
    RCP<MulitiVector> owned_solution;
    RCP<Vector> owned_residual;
    RCP<Matrix> owned_jacobian;
    RCP<MultiVector> ovlp_solution;
    RCP<Vector> ovlp_residual;
    RCP<Vector> ovlp_jacobain;
};

}

#endif
