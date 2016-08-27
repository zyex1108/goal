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
    unsigned num_eqs;

    std::vector<PHX::MDField<ScalarT, Elem, Node> > z;

PHX_EVALUATOR_CLASS_END

}

#endif
