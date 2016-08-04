#ifndef goal_ev_qoi_avg_displacement_hpp
#define goal_ev_qoi_avg_displacement_hpp

#include "phx_macros.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

struct Layouts;

PHX_EVALUATOR_CLASS(QoIAvgDisplacement)

  private:

    Teuchos::RCP<Layouts> dl;
    Teuchos::RCP<const ParameterList> params;

    unsigned num_dims;
    unsigned num_nodes;
    unsigned num_qps;

    Teuchos::Array<std::string> disp_names;

    PHX::MDField<double, Elem, QP> wDv;
    std::vector<PHX::MDField<ScalarT, Elem, QP> > disp;
    PHX::MDField<ScalarT, Elem> avg_disp;

PHX_EVALUATOR_CLASS_END

}

#endif
