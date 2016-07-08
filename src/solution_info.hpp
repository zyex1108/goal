#ifndef goal_solution_info_hpp
#define goal_solution_info_hpp

#include "data_types.hpp"
#include "Teuchos_RCP.hpp"

namespace goal {

using Teuchos::rcp;
using Teuchos::RCP;

class Mesh;

class SolutionInfo
{
  public:
    void resize(RCP<Mesh> m, bool enable_dynamics);
    void gather_solution();
    void scatter_solution();
    void gather_residual();
    void gather_jacobian();
    RCP<MultiVector> owned_solution;
    RCP<Vector> owned_residual;
    RCP<Matrix> owned_jacobian;
    RCP<MultiVector> ovlp_solution;
    RCP<Vector> ovlp_residual;
    RCP<Matrix> ovlp_jacobian;
    RCP<Export> exporter;
    RCP<Import> importer;
};

RCP<SolutionInfo> sol_info_create(RCP<Mesh> m, bool enable_dynamics);

}

#endif
