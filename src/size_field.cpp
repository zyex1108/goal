#include "size_field.hpp"
#include "mesh.hpp"
#include "assert_param.hpp"

#include <ma.h>

namespace goal {

void validate_size_params(
    RCP<const ParameterList> p,
    RCP<Mesh> m)
{
}

ma::IsotropicFunction* get_size_field(
    RCP<const ParameterList> p,
    RCP<Mesh> m)
{
  return 0;
}

}
