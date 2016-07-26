#ifndef goal_ev_elem_size_hpp
#define goal_ev_elem_size_hpp

#include "phx_macros.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

struct Layouts;

PHX_EVALUATOR_CLASS(ElemSize)

  private:

    RCP<Layouts> dl;

    unsigned num_nodes;
    unsigned num_qps;
    unsigned num_dims;

    PHX::MDField<double, Elem, Node, QP, Dim> gBF;
    PHX::MDField<double, Elem, QP> size;

PHX_EVALUATOR_CLASS_END

}

#endif
