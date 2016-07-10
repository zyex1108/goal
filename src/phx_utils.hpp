#ifndef goal_phx_utils_hpp
#define goal_phx_utils_hpp

#include "dimension.hpp"
#include "layouts.hpp"
#include <Phalanx.hpp>

namespace goal {

template <typename ScalarT>
void get_field(
    std::string const& name,
    RCP<Layouts> dl,
    PHX::MDField<ScalarT,Elem,Node>& f)
{
  f = PHX::MDField<ScalarT,Elem,Node>(name,dl->node_scalar);
}

template <typename ScalarT>
void get_grad_field(
    PHX::MDField<ScalarT, Elem, Node, QP> f,
    RCP<Layouts> dl,
    PHX::MDField<ScalarT, Elem, Node, QP, Dim>& gf)
{
  std::string const& name = f.fieldTag().name();
  gf = PHX::MDField<ScalarT, Elem, Node, QP, Dim>
    ("grad_" + name, dl->node_qp_vector);
}

}

#endif
