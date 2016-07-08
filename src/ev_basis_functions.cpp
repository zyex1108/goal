#include "ev_basis_functions.hpp"
#include "mesh.hpp"
#include "layouts.hpp"
#include "workset.hpp"

#include <apf.h>
#include <apfMesh2.h>

namespace goal {


PHX_EVALUATOR_CTOR(BasisFunctions, p)
{
}

PHX_POST_REGISTRATION_SETUP(BasisFunctions, data, fm)
{
}

PHX_EVALUATE_FIELDS(BasisFunctions, workset)
{
}

}
