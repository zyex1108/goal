#ifndef goal_size_field_hpp
#define goal_size_field_hpp

#include <Teuchos_RCP.hpp>

namespace Teuchos {
class ParameterList;
}

namespace ma {
class IsotropicFunction;
}

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class Mesh;

void validate_size_params(
    RCP<const ParameterList> p,
    RCP<Mesh> mesh);

ma::IsotropicFunction* get_size_field(
    RCP<const ParameterList> p,
    RCP<Mesh> mesh);

}

#endif
