#ifndef goal_workset_hpp
#define goal_workset_hpp

#include "data_types.hpp"
#include <Teuchos_RCP.hpp>

namespace apf {
class MeshEntity;
}

namespace goal {

using Teuchos::RCP;

struct Workset
{
  unsigned size;
  std::string set;
  double t_new;
  double t_old;
  RCP<MultiVector> u;
  RCP<Vector> r;
  RCP<Matrix> J;
  double alpha;
  double beta;
  double gamma;
  std::vector<apf::MeshEntity*> ents;
  bool is_adjoint;
  RCP<Vector> z;
  RCP<Vector> q;
};

}

#endif
