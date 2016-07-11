#ifndef goal_ev_kinematics_hpp
#define goal_ev_kinematics_hpp

#include "phx_macros.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

struct Layouts;

PHX_EVALUATOR_CLASS(Kinematics)

  private:

    RCP<Layouts> dl;
    Teuchos::Array<std::string> disp_names;

    unsigned num_qps;
    unsigned num_dims;

    std::vector<PHX::MDField<ScalarT, Elem, QP, Dim> > grad_u;

    PHX::MDField<ScalarT, Elem, QP, Dim, Dim> def_grad;
    PHX::MDField<ScalarT, Elem, QP> det_def_grad;

PHX_EVALUATOR_CLASS_END

}

#endif
