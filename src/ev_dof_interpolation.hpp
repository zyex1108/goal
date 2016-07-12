#ifndef goal_dof_interpolation_hpp
#define goal_dof_interpolation_hpp

#include "phx_macros.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

struct Layouts;

PHX_EVALUATOR_CLASS(DOFInterpolation)

  private:

    RCP<Layouts> dl;
    Teuchos::Array<std::string> names;
    unsigned index;

    unsigned num_nodes;
    unsigned num_qps;
    unsigned num_dims;
    unsigned num_eqs;

    PHX::MDField<double, Elem, Node, QP> BF;
    PHX::MDField<double, Elem, Node, QP, Dim> gBF;
    std::vector<PHX::MDField<ScalarT, Elem, Node> > nodal;

    std::vector<PHX::MDField<ScalarT, Elem, QP> > dofs;
    std::vector<PHX::MDField<ScalarT, Elem, QP, Dim> > gdofs;

PHX_EVALUATOR_CLASS_END

}

#endif
