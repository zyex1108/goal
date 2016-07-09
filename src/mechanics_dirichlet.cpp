#include "mechanics.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "mesh.hpp"

namespace goal {

template <typename EvalT>
void goal::Mechanics::register_dirichlet(FieldManager fm)
{
}

}

template void goal::Mechanics::
register_dirichlet<goal::GoalTraits::Residual>(FieldManager fm);

template void goal::Mechanics::
register_dirichlet<goal::GoalTraits::Jacobian>(FieldManager fm);
