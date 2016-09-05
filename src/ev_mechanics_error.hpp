#ifndef goal_mechanics_error_hpp
#define goal_mechanics_error_hpp

#include "phx_macros.hpp"
#include "traits.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::ParameterList;

class Mesh;
struct Layouts;

PHX_EVALUATOR_CLASS(MechanicsError)

  private:

    Teuchos::RCP<Layouts> dl;
    Teuchos::RCP<Mesh> mesh;

    unsigned num_nodes;
    unsigned num_qps;
    unsigned num_dims;

    PHX::MDField<ScalarT, Elem, QP, Dim, Dim> first_pk;

PHX_EVALUATOR_CLASS_END

}

#endif
