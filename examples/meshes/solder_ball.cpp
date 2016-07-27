#include <gmodel.hpp>

using namespace gmod;

static double square(double x) {return x * x;}

static ObjPtr new_solder_ball(
    Vector center,
    double radius,
    double thickness)
{
  PointPtr cen_pt = new_point2(center);
  ObjPtr shell = new_shell();
  double c = sqrt(square(radius) - square(thickness / 2));
  double d = 1.0;
  ObjPtr circles[2];
  for (unsigned i = 0; i < 2; ++i) {
    Vector circ_cen = add_vectors(center, Vector{0, 0, d * thickness / 2});
    circles[i] = new_circle(circ_cen, Vector{0,0,1}, Vector{c,0,0});
    add_use(shell, i, new_plane2(circles[i]));
    d = -d;
  }
  std::vector<PointPtr> circle_points[2] = {
    loop_points(circles[0]),
    loop_points(circles[1])};
  ObjPtr verticals[4];
  for (unsigned i = 0; i < 4; ++i) {
    verticals[i] = new_arc2(
        circle_points[1][i],
        cen_pt,
        circle_points[0][i]);
  }
  for (unsigned i = 0; i < 4; ++i) {
    ObjPtr loop = new_loop();
    add_use(loop, FORWARD, circles[1]->used[i].obj);
    add_use(loop, FORWARD, verticals[(i + 1) % 4]);
    add_use(loop, REVERSE, circles[0]->used[i].obj);
    add_use(loop, REVERSE, verticals[i]);
    add_use(shell, FORWARD, new_ruled2(loop));
  }
  return new_volume2(shell);
}

static ObjPtr ball_bottom(ObjPtr ball)
{
  return volume_shell(ball)->used[1].obj;
}

static ObjPtr ball_top(ObjPtr ball)
{
  return volume_shell(ball)->used[0].obj;
}

static double const top_thickness = 1.0;
static double const bottom_thickness = 1.0;
static double const ball_thickness = 0.5;
static double const ball_radius = 0.4;
static double const unit_width = 1.0;
static unsigned const nunits_per_side = 2;

int main()
{
  default_size = 0.3;
  ObjPtr g = new_group();
  ObjPtr bot = new_cube(
      Vector{0,0,0},
      Vector{unit_width * nunits_per_side, 0, 0},
      Vector{0, unit_width * nunits_per_side, 0},
      Vector{0, 0, bottom_thickness});
  add_to_group(g, bot);
  ObjPtr top = new_cube(
      Vector{0, 0, bottom_thickness + ball_thickness},
      Vector{unit_width * nunits_per_side, 0, 0},
      Vector{0, unit_width * nunits_per_side, 0},
      Vector{0, 0, top_thickness});
  add_to_group(g, top);
  default_size = 0.1;
  for (unsigned i = 0; i < nunits_per_side; ++i)
    for (unsigned j = 0; j < nunits_per_side; ++j) {
      ObjPtr sb = new_solder_ball(
          Vector{
          i * unit_width + (unit_width / 2),
          j * unit_width + (unit_width / 2),
          bottom_thickness + (ball_thickness / 2)},
          ball_radius, ball_thickness);
      add_to_group(g, sb);
      weld_volume_face_into(bot, sb, get_cube_face(bot, TOP), ball_bottom(sb));
      weld_volume_face_into(top, sb, get_cube_face(top, BOTTOM), ball_top(sb));
    }
  write_closure_to_geo(g, "solder_ball.geo");
  write_closure_to_dmg(g, "solder_ball.dmg");
}
