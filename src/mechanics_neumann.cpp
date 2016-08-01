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

/* ETI */

template void goal::Mechanics::
register_neumann<goal::GoalTraits::Forward>(FieldManager fm);

template void goal::Mechanics::
register_neumann<goal::GoalTraits::Derivative>(FieldManager fm);
