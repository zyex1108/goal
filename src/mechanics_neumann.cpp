#include "mechanics.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "mesh.hpp"

namespace goal {

template <typename EvalT>
void goal::Mechanics::register_neumann(FieldManager fm)
{
}

}

template void goal::Mechanics::
register_neumann<goal::GoalTraits::Residual>(FieldManager fm);

template void goal::Mechanics::
register_neumann<goal::GoalTraits::Jacobian>(FieldManager fm);
