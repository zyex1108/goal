#include "assert_param.hpp"
#include "control.hpp"

#include <Teuchos_ParameterList.hpp>

namespace goal {

void assert_param(RCP<const ParameterList> p, std::string const& v)
{
  if (! p->isParameter(v))
    fail("the parameter '%s' was expected to be a member "
        "of the parameterlist '%s', but it was not found",
        v.c_str(), p->name().c_str());
}

void assert_sublist(RCP<const ParameterList> p, std::string const& v)
{
  if (! p->isSublist(v))
    fail("the sublist '%s' was expected to be a member "
        "of the parameterlist '%s', but it was not found",
        v.c_str(), p->name().c_str());
}

}
