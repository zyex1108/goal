#include "output.hpp"
#include "mesh.hpp"
#include "mechanics.hpp"
#include "solution_info.hpp"
#include "control.hpp"
#include "data_types.hpp"

#include <PCU.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>

namespace goal {

static RCP<ParameterList> get_valid_params()
{
  RCP<ParameterList> p = rcp(new ParameterList);
  p->set<std::string>("out file", "");
  p->set<bool>("turn off", false);
  p->set<unsigned>("output interval", 0);
  p->set<bool>("save stabilized stress", false);
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
  index(0),
  save_stabilized(false)
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
  if (params->isParameter("save stabilized stress"))
    save_stabilized = params->get<bool>("save stabilized stress");
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

static void attach_vector(
    RCP<const Vector> u,
    RCP<Mesh> mesh,
    RCP<Mechanics> mech,
    Teuchos::Array<std::string> const& names)
{
  ArrayRCP<const ST> data = u->get1dView();
  apf::Mesh* m = mesh->get_apf_mesh();
  apf::DynamicArray<apf::Node> nodes = mesh->get_apf_nodes();
  std::vector<apf::Field*> fields;
  for (unsigned j=0; j < mesh->get_num_eqs(); ++j) {
    apf::Field* f = apf::createFieldOn(m, names[j].c_str(), apf::SCALAR);
    fields.push_back(f);
  }
  for (unsigned i=0; i < nodes.size(); ++i) {
    apf::Node* node = &(nodes[i]);
    if (! m->isOwned(node->entity)) continue;
    if (m->getType(node->entity) != apf::Mesh::VERTEX) continue;
    for (unsigned j=0; j < names.size(); ++j) {
      unsigned eq = mech->get_offset(names[j]);
      LO row = mesh->get_lid(node, eq);
      double v = data[row];
      apf::setScalar(fields[j], node->entity, node->node, v);
    }
  }
  for (unsigned j=0; j < names.size(); ++j)
    apf::synchronize(fields[j]);
}

static void attach_solutions(
    RCP<Mesh> mesh,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s)
{
  RCP<MultiVector> sv = s->owned_solution;
  unsigned nv = sv->getNumVectors();
  for (unsigned i=0; i < nv; ++i) {
    RCP<const Vector> u = sv->getVector(i);
    Teuchos::Array<std::string> names = mech->get_var_names(i);
    attach_vector(u, mesh, mech, names);
  }
}

static void destroy_solutions(
    RCP<Mesh> mesh,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s)
{
  apf::Mesh* m = mesh->get_apf_mesh();
  RCP<MultiVector> sv = s->owned_solution;
  unsigned nv = sv->getNumVectors();
  for (unsigned i=0; i < nv; ++i) {
    Teuchos::Array<std::string> names = mech->get_var_names(i);
    for (unsigned j=0; j < names.size(); ++j) {
      apf::Field* f = m->findField(names[j].c_str());
      CHECK(f);
      apf::destroyField(f);
    }
  }
}

static void stabilize(RCP<Mesh> mesh)
{
  apf::Mesh2* m = mesh->get_apf_mesh();
  apf::Field* pressure = m->findField("p");
  apf::Field* cauchy = m->findField("cauchy");
  CHECK(pressure);
  CHECK(cauchy);
  unsigned num_dims = mesh->get_num_dims();
  unsigned num_qps = mesh->get_num_elem_qps();
  unsigned order = mesh->get_q_order();
  apf::Vector3 point;
  apf::Matrix3x3 sigma;
  apf::MeshEntity* elem;
  apf::MeshIterator* elems = m->begin(m->getDimension());
  while ((elem = m->iterate(elems))) {
    apf::MeshElement* me = apf::createMeshElement(m, elem);
    apf::Element* pe = apf::createElement(pressure, me);
    for (unsigned qp=0; qp < num_qps; ++qp) {
      apf::getIntPoint(me, order, qp, point);
      double p = apf::getScalar(pe, point);
      apf::getMatrix(cauchy, elem, qp, sigma);
      double dilatation = sigma[0][0];
      for (unsigned i=1; i < num_dims; ++i)
        dilatation += sigma[i][i];
      dilatation /= num_dims;
      for (unsigned i=0; i < num_dims; ++i)
        sigma[i][i] += p - dilatation;
      apf::setMatrix(cauchy, elem, qp, sigma);
    }
    apf::destroyElement(pe);
    apf::destroyMeshElement(me);
  }
  m->end(elems);
}

void Output::write(const double t)
{
  if (turn_off) return;
  static unsigned my_out_interval = 0;
  if (my_out_interval++ % out_interval) return;
  attach_solutions(mesh, mechanics, sol_info);
  if (save_stabilized) stabilize(mesh);
  write_vtk(t);
  destroy_solutions(mesh, mechanics, sol_info);
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
