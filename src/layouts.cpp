#include "layouts.hpp"
#include "dimension.hpp"

namespace goal {

using PHX::MDALayout;

Layouts::Layouts()
{
  elem_scalar = rcp(new MDALayout<Dummy>(0));
  node_scalar = rcp(new MDALayout<Dummy>(0));
  node_vector = rcp(new MDALayout<Dummy>(0));
  qp_scalar = rcp(new MDALayout<Dummy>(0));
  qp_vector = rcp(new MDALayout<Dummy>(0));
  qp_tensor = rcp(new MDALayout<Dummy>(0));
  node_qp_scalar = rcp(new MDALayout<Dummy>(0));
  node_qp_vector = rcp(new MDALayout<Dummy>(0));
  dummy = rcp(new MDALayout<Dummy>(0));
}

Layouts::Layouts(size_t ws_size, size_t n_nodes, size_t n_qps, size_t n_dims)
{
  elem_scalar = rcp(new MDALayout<Elem>(ws_size));
  node_scalar = rcp(new MDALayout<Elem, Node>(ws_size, n_nodes));
  node_vector = rcp(new MDALayout<Elem, Node, Dim>(ws_size, n_nodes, n_dims));
  qp_scalar = rcp(new MDALayout<Elem, QP>(ws_size, n_qps));
  qp_vector = rcp(new MDALayout<Elem, QP, Dim>(ws_size, n_qps, n_dims));
  qp_tensor = rcp(new MDALayout<Elem, QP, Dim, Dim>(ws_size, n_qps, n_dims, n_dims));
  node_qp_scalar = rcp(new MDALayout<Elem, Node, QP>(ws_size, n_nodes, n_qps));
  node_qp_vector = rcp(new MDALayout<Elem, Node, QP, Dim>(ws_size, n_nodes, n_qps, n_dims));
  dummy = rcp(new MDALayout<Dummy>(0));
}

}
