#ifndef goal_ev_mechanics_residual_hpp
#define goal_ev_mechanics_residual_hpp

#include "phx_macros.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

struct Layouts;
class Mesh;

PHX_EVALUATOR_CLASS(MechanicsResidual)

  private:

    RCP<Layouts> dl;
    RCP<Mesh> mesh;

    std::string bf_name;
    Teuchos::Array<std::string> disp_names;

    bool enable_dynamics;
    bool have_body_force;

    double rho;
    Teuchos::Array<std::string> body_force;

    unsigned num_nodes;
    unsigned num_qps;
    unsigned num_dims;

    PHX::MDField<double, Elem, QP> wDv;
    PHX::MDField<double, Elem, Node, QP> BF;
    PHX::MDField<double, Elem, Node, QP, Dim> gBF;
    PHX::MDField<ScalarT, Elem, QP, Dim, Dim> stress;
    std::vector<PHX::MDField<ScalarT, Elem, QP> > acc;
    std::vector<PHX::MDField<ScalarT, Elem, Node> > resid;

PHX_EVALUATOR_CLASS_END

}

#endif
