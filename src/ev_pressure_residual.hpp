#ifndef goal_ev_pressure_residual_hpp
#define goal_ev_pressure_residual_hpp

#include "phx_macros.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

struct Layouts;

PHX_EVALUATOR_CLASS(PressureResidual)

  private:

    RCP<Layouts> dl;
    std::string bf_name;
    std::string pressure_name;

    unsigned num_nodes;
    unsigned num_qps;
    unsigned num_dims;

    PHX::MDField<double, Elem, QP> wDv;
    PHX::MDField<double, Elem, QP> BF;
    PHX::MDField<double, Elem, QP, Dim> gBF;
    PHX::MDField<ScalarT, Elem, QP, Dim, Dim> def_grad;
    PHX::MDField<ScalarT, Elem, QP, Dim, Dim> stress;
    PHX::MDField<ScalarT, Elem, QP> pressure;
    PHX::MDField<ScalarT, Elem, QP, Dim> pressure_grad;
    PHX::MDField<ScalarT, Elem, Node> resid;

PHX_EVALUATOR_CLASS_END

}

#endif
