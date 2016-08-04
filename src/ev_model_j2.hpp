#ifndef goal_ev_model_j2_hpp
#define goal_ev_model_j2_hpp

#include "phx_macros.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

struct Layouts;
class StateFields;

PHX_EVALUATOR_CLASS(ModelJ2)

  private:

    RCP<Layouts> dl;
    RCP<StateFields> states;
    RCP<const ParameterList> params;

    double E; /* elastic modulus */
    double nu; /* poisson's ratio */
    double K; /* hardening modulus */
    double Y; /* yield strength */

    unsigned num_nodes;
    unsigned num_qps;
    unsigned num_dims;

    PHX::MDField<ScalarT, Elem, QP, Dim, Dim> def_grad;
    PHX::MDField<ScalarT, Elem, QP> det_def_grad;
    PHX::MDField<ScalarT, Elem, QP, Dim, Dim> stress;

PHX_EVALUATOR_CLASS_END

}

#endif
