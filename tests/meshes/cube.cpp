#include "gmodel.hpp"

int main()
{
  gmod::default_size = 0.9;
  auto s = gmod::new_cube(
      gmod::Vector{0,0,0},
      gmod::Vector{1,0,0},
      gmod::Vector{0,1,0},
      gmod::Vector{0,0,1});
  gmod::write_closure_to_geo(s, "cube.geo");
  gmod::write_closure_to_dmg(s, "cube.dmg");
}
