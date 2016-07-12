#ifndef goal_layouts_hpp
#define goal_layouts_hpp

#include <Teuchos_RCP.hpp>
#include <Phalanx_DataLayout_MDALayout.hpp>

namespace goal {

using Teuchos::rcp;
using Teuchos::RCP;

struct Layouts
{
  Layouts();
  Layouts(size_t ws_size, size_t n_nodes, size_t n_qps, size_t n_dims);
  RCP<PHX::DataLayout> elem_scalar;
  RCP<PHX::DataLayout> node_scalar;
  RCP<PHX::DataLayout> node_vector;
  RCP<PHX::DataLayout> qp_scalar;
  RCP<PHX::DataLayout> qp_vector;
  RCP<PHX::DataLayout> qp_tensor;
  RCP<PHX::DataLayout> node_qp_scalar;
  RCP<PHX::DataLayout> node_qp_vector;
  RCP<PHX::DataLayout> dummy;
};

}

#endif
