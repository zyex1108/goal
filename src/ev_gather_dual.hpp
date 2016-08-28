#ifndef goal_gather_dual_hpp
#define goal_gather_dual_hpp

#include "phx_macros.hpp"
#include "traits.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::ParameterList;

class Mesh;
class Layouts;

PHX_EVALUATOR_CLASS(GatherDual)

  private:

    Teuchos::RCP<Layouts> dl;
    Teuchos::RCP<Mesh> mesh;
    Teuchos::Array<std::string> names;

    unsigned num_nodes;
    unsigned num_qps;
    unsigned num_dims;
    unsigned num_eqs;

    PHX::MDField<double, Elem, Node, QP> BF;
    PHX::MDField<double, Elem, Node, QP, Dim> gBF;
    std::vector<PHX::MDField<ScalarT, Elem, Node> > nodal;
    std::vector<PHX::MDField<ScalarT, Elem, QP> > duals;
    std::vector<PHX::MDField<ScalarT, Elem, QP, Dim> > gduals;

PHX_EVALUATOR_CLASS_END

}

#endif
