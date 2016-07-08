#ifndef goal_assert_param_hpp
#define goal_assert_param_hpp

#include <Teuchos_RCP.hpp>

namespace Teuchos {
class ParameterList;
}

namespace goal {

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::ParameterList;

void assert_param(RCP<const ParameterList> p, std::string const& v);

}

#endif
