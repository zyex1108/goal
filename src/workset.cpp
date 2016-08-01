#include "workset.hpp"

namespace goal {

Workset::Workset() :
  size(0),
  set(""),
  t_new(0.0),
  t_old(0.0),
  alpha(0.0),
  beta(0.0),
  gamma(0.0),
  is_adjoint(false)
{
}

}
