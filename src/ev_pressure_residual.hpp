#ifndef goal_ev_pressure_residual_hpp
#define goal_ev_pressure_residual_hpp

#include "phx_macros.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

struct Layouts;
class Mesh;

PHX_EVALUATOR_CLASS(PressureResidual)

  private:

    RCP<Layouts> dl;
    RCP<const ParameterList> params;
    bool small_strain;
    std::string pressure_name;

    double G;
    double K;
    double alpha;

    unsigned num_nodes;
    unsigned num_qps;
    unsigned num_dims;

    PHX::MDField<double, Elem, QP> size;
    PHX::MDField<double, Elem, QP> wDv;
    PHX::MDField<double, Elem, Node, QP> BF;
    PHX::MDField<double, Elem, Node, QP, Dim> gBF;
    PHX::MDField<ScalarT, Elem, QP, Dim, Dim> def_grad;
    PHX::MDField<ScalarT, Elem, QP> det_def_grad;
    PHX::MDField<ScalarT, Elem, QP, Dim, Dim> stress;
    PHX::MDField<ScalarT, Elem, QP> pressure;
    PHX::MDField<ScalarT, Elem, QP, Dim> pressure_grad;
    PHX::MDField<ScalarT, Elem, Node> resid;

PHX_EVALUATOR_CLASS_END

}

#endif
