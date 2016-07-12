#ifndef goal_ev_basis_functions_hpp
#define goal_ev_basis_functions_hpp

#include "phx_macros.hpp"

namespace goal {

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::ParameterList;

class Mesh;
struct Layouts;

PHX_EVALUATOR_CLASS(BasisFunctions)

  private:

    RCP<Layouts> dl;
    RCP<Mesh> mesh;

    unsigned num_nodes;
    unsigned num_qps;
    unsigned num_dims;

    unsigned p_order;
    unsigned q_order;

    PHX::MDField<double, Elem, QP> wDv;
    PHX::MDField<double, Elem, Node, QP> BF;
    PHX::MDField<double, Elem, Node, QP, Dim> gBF;

PHX_EVALUATOR_CLASS_END

}

#endif
