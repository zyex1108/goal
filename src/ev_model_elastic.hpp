#ifndef goal_ev_model_elastic_hpp
#define goal_ev_model_elastic_hpp

#include "phx_macros.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

struct Layouts;
class StateFields;

PHX_EVALUATOR_CLASS(ModelElastic)

  private:

    RCP<Layouts> dl;
    RCP<StateFields> states;
    RCP<const ParameterList> params;
    RCP<const ParameterList> temp_params;
    Teuchos::Array<std::string> disp_names;

    bool have_temp;

    double E; /* elastic modulus */
    double nu; /* poisson's ratio */
    double alpha; /* thermal expansion coefficient */

    unsigned num_nodes;
    unsigned num_qps;
    unsigned num_dims;

    std::vector<PHX::MDField<ScalarT, Elem, QP, Dim> > grad_u;
    PHX::MDField<ScalarT, Elem, QP, Dim, Dim> stress;

PHX_EVALUATOR_CLASS_END

}

#endif
