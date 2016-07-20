#ifndef goal_ev_first_pk_hpp
#define goal_ev_first_pk_hpp

#include "phx_macros.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

struct Layouts;

PHX_EVALUATOR_CLASS(FirstPK)

  private:

    RCP<Layouts> dl;

    unsigned num_qps;
    unsigned num_dims;

    bool small_strain;
    bool have_pressure;

    PHX::MDField<ScalarT, Elem, QP, Dim, Dim> def_grad;
    PHX::MDField<ScalarT, Elem, QP> det_def_grad;
    PHX::MDField<ScalarT, Elem, QP> pressure;
    PHX::MDField<ScalarT, Elem, QP, Dim, Dim> cauchy;

    PHX::MDField<ScalarT, Elem, QP, Dim, Dim> first_pk;

PHX_EVALUATOR_CLASS_END

}

#endif
