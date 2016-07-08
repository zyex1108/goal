#ifndef troy_workset_hpp
#define troy_workset_hpp

#include "data_types.hpp"
#include <Teuchos_RCP.hpp>

namespace apf {
class MeshEntity;
}

namespace goal {

using Teuchos::RCP;

struct workset
{
  unsigned size;
  std::string set;
  double t_new;
  double t_old;
  RCP<MultiVector> u;
  RCP<Vector> r;
  RCP<Matrix> J;
  std::vector<apf::MeshEntity*> ents;
};

}

#endif
