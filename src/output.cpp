#include "output.hpp"
#include "mesh.hpp"
#include "mechanics.hpp"
#include "control.hpp"
#include "data_types.hpp"

#include <PCU.h>
#include <apf.h>
#include <apfMesh2.h>

namespace goal {

static RCP<ParameterList> get_valid_params()
{
  RCP<ParameterList> p = rcp(new ParameterList);
  p->set<std::string>("out file", "");
  p->set<bool>("turn off", false);
  p->set<unsigned>("output interval", 0);
  return p;
}

static void validate_params(RCP<const ParameterList> p)
{
  CHECK(p->isParameter("out file"));
  p->validateParameters(*get_valid_params(), 0);
}

static void write_initial_pvd(
    std::string const& name,
    unsigned& o_file_pos)
{
  if (! PCU_Comm_Self()) {
    std::string pvd = name + ".pvd";
    std::fstream pvdf;
    pvdf.open(pvd.c_str(), std::ios::out);
    pvdf << "<VTKFile type=\"Collection\" version=\"0.1\">" << std::endl;
    pvdf << "  <Collection>" << std::endl;
    o_file_pos = pvdf.tellp();
    pvdf << "  </Collection>" << std::endl;
    pvdf << "</VTKFile>" << std::endl;
    pvdf.close();
  }
}

Output::Output(
    RCP<const ParameterList> p,
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s) :
  params(p),
  mesh(m),
  mechanics(mech),
  sol_info(s),
  turn_off(false),
  out_interval(1),
  o_file_pos(0),
  index(0)
{
  validate_params(params);
  get_inputs();
  if (! turn_off)
    write_initial_pvd(name, o_file_pos);
}

void Output::get_inputs()
{
  name = params->get<std::string>("out file");
  if (params->isParameter("turn off"))
    turn_off = params->get<bool>("turn off");
  if (params->isParameter("output interval"))
    out_interval = params->get<unsigned>("output interval");
}

static void update_pvd(
    std::string const& name,
    std::string const& vtu,
    unsigned& o_file_pos,
    const double t)
{
  if (! PCU_Comm_Self()) {
    std::string pvd = name + ".pvd";
    std::fstream pvdf;
    pvdf.open(pvd.c_str(), std::ios::out | std::ios::in);
    pvdf.seekp(o_file_pos);
    pvdf << "    <DataSet timestep=\"" << t << "\" group=\"\" ";
    pvdf << "part=\"0\" file=\"" << vtu << "/" << vtu;
    pvdf << ".pvtu\"/>" << std::endl;
    o_file_pos = pvdf.tellp();
    pvdf << "  </Collection>" << std::endl;
    pvdf << "</VTKFile>" << std::endl;
    pvdf.close();
  }
}

void Output::write_vtk(const double t)
{
  std::ostringstream oss;
  oss << name << "_" << index;
  std::string vtu = oss.str();
  update_pvd(name, vtu, o_file_pos, t);
  apf::writeVtkFiles(vtu.c_str(), mesh->get_apf_mesh());
  ++index;
}

void Output::write(const double t)
{
  if (turn_off) return;
  static unsigned my_out_interval = 0;
  if (my_out_interval++ % out_interval) return;
  write_vtk(t);
}

RCP<Output> output_create(
    RCP<const ParameterList> p,
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s)
{
  RCP<const ParameterList> op = rcpFromRef(p->sublist("output"));
  return rcp(new Output(op, m, mech, s));
}

}
