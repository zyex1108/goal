#include "gmodel.hpp"

int main()
{
  auto s = gmod::new_square(
      gmod::Vector{0,0,0},
      gmod::Vector{1,0,0},
      gmod::Vector{0,1,0});
  gmod::write_closure_to_geo(s, "square.geo");
  gmod::write_closure_to_dmg(s, "square.dmg");
}
