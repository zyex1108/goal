#include "gmodel.hpp"

int main()
{
  gmod::default_size = 0.15;
  auto s = gmod::new_square(
      gmod::Vector{0,0,0},
      gmod::Vector{5,0,0},
      gmod::Vector{0,1,0});
  gmod::write_closure_to_geo(s, "rectangle.geo");
  gmod::write_closure_to_dmg(s, "rectangle.dmg");
}
