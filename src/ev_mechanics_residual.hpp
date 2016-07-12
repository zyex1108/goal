#ifndef goal_ev_mechanics_residual_hpp
#define goal_ev_mechanics_residual_hpp

#include "phx_macros.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

struct Layouts;

PHX_EVALUATOR_CLASS(MechanicsResidual)

  private:

    RCP<Layouts> dl;
    std::string bf_name;
    Teuchos::Array<std::string> disp_names;

    unsigned num_nodes;
    unsigned num_qps;
    unsigned num_dims;

    PHX::MDField<double, Elem, QP> wDv;
    PHX::MDField<double, Elem, Node, QP, Dim> gBF;
    PHX::MDField<ScalarT, Elem, QP, Dim, Dim> stress;
    std::vector<PHX::MDField<ScalarT, Elem, Node> > resid;

PHX_EVALUATOR_CLASS_END

}

#endif
