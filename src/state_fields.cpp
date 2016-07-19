#include "state_fields.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>

namespace goal {

const apf::ValueType apf_types[3] =
{ apf::SCALAR,
  apf::VECTOR,
  apf::MATRIX };

StateFields::StateFields(RCP<Mesh> m) :
  mesh(m),
  apf_mesh(m->get_apf_mesh())
{
}

static void set_identity(RCP<Mesh> m, apf::Field* f, FieldType t)
{
  CHECK(t == TENSOR);
  apf::Matrix3x3 I(1,0,0,0,1,0,0,0,1);
  apf::Mesh* am = m->get_apf_mesh();
  apf::MeshEntity* elem;
  apf::MeshIterator* it = am->begin(m->get_num_dims());
  while ((elem = am->iterate(it)))
    for (unsigned qp=0; qp < m->get_num_elem_qps(); ++qp)
      apf::setMatrix(f, elem, qp, I);
  am->end(it);
}

void StateFields::add(const char* n, FieldType t, bool save, bool I)
{
  unsigned dim = mesh->get_num_dims();
  unsigned q_order = mesh->get_q_order();
  apf::FieldShape* s = apf::getIPFitShape(dim, q_order);
  apf::Field* f = apf::createField(apf_mesh, n, apf_types[t], s);
  if (! I) apf::zeroField(f);
  else set_identity(mesh, f, t);
  states.push_back(f);
  if (save) {
    std::string oname = std::string(n) + "_old";
    apf::Field* g = apf::createField(
        apf_mesh, oname.c_str(), apf_types[t], s);
    if (!I) apf::zeroField(g);
    else set_identity(mesh, g, t);
    states.push_back(g);
    old_states.push_back(g);
  }
}

static double get_val(double v)
{
  return v;
}

static double get_val(FadType const& v)
{
  return v.val();
}

static void zero(apf::Vector3& v)
{
  for (unsigned i=0; i < 3; ++i)
    v[i] = 0.0;
}

static void zero(apf::Matrix3x3& m)
{
  for (unsigned i=0; i < 3; ++i)
  for (unsigned j=0; j < 3; ++j)
    m[i][j] = 0.0;
}

template <typename T>
void StateFields::set_scalar(
    char const* name,
    apf::MeshEntity* e,
    unsigned n,
    T const& v)
{
  apf::Field* f = apf_mesh->findField(name);
  apf::setScalar(f, e, n, get_val(v));
}

template <typename T>
void StateFields::set_vector(
    char const* name,
    apf::MeshEntity* e,
    unsigned n,
    Intrepid2::Vector<T> const& v)
{
  apf::Vector3 val;
  zero(val);
  for (unsigned i=0; i < v.get_dimension(); ++i)
    val[i] = get_val(v(i));
  apf::Field* f = apf_mesh->findField(name);
  apf::setVector(f, e, n, val);
}

template <typename T>
void StateFields::set_tensor(
    char const* name,
    apf::MeshEntity* e,
    unsigned n,
    Intrepid2::Tensor<T> const& v)
{
  apf::Matrix3x3 val;
  zero(val);
  for (unsigned i=0; i < v.get_dimension(); ++i)
  for (unsigned j=0; j < v.get_dimension(); ++j)
    val[i][j] = get_val(v(i,j));
  apf::Field* f = apf_mesh->findField(name);
  apf::setMatrix(f, e, n, val);

}

template <typename T>
void StateFields::get_scalar(
    char const* name,
    apf::MeshEntity* e,
    unsigned n,
    T& v)
{
  apf::Field* f = apf_mesh->findField(name);
  double val = apf::getScalar(f, e, n);
  v = (T)val;
}

template <typename T>
void StateFields::get_vector(
    char const* name,
    apf::MeshEntity* e,
    unsigned n,
    Intrepid2::Vector<T>& v)
{
  apf::Field* f = apf_mesh->findField(name);
  apf::Vector3 val;
  apf::getVector(f, e, n, val);
  for (unsigned i=0; i < v.get_dimension(); ++i)
    v(i) = (T)val[i];
}

template <typename T>
void StateFields::get_tensor(
    char const* name,
    apf::MeshEntity* e,
    unsigned n,
    Intrepid2::Tensor<T>& v)
{
  apf::Field* f = apf_mesh->findField(name);
  apf::Matrix3x3 val;
  apf::getMatrix(f, e, n, val);
  for (unsigned i=0; i < v.get_dimension(); ++i)
  for (unsigned j=0; j < v.get_dimension(); ++j)
    v(i,j) = (T)val[i][j];
}

void StateFields::update()
{
  for (unsigned i=0; i < old_states.size(); ++i) {
    apf::Field* to = old_states[i];
    std::string to_name = (std::string)apf::getName(to);
    std::string from_name = to_name.erase(to_name.find("_old"),4);
    apf::Field* from = apf_mesh->findField(from_name.c_str());
    apf::copyData(to,from);
  }
}

/* ETI */
template void StateFields::set_scalar(char const* name, apf::MeshEntity* e, unsigned n, double const& v);
template void StateFields::set_scalar(char const* name, apf::MeshEntity* e, unsigned n, FadType const& v);
template void StateFields::set_vector(char const* name, apf::MeshEntity* e, unsigned n, Intrepid2::Vector<double> const& v);
template void StateFields::set_vector(char const* name, apf::MeshEntity* e, unsigned n, Intrepid2::Vector<FadType> const& v);
template void StateFields::set_tensor(char const* name, apf::MeshEntity* e, unsigned n, Intrepid2::Tensor<double> const& v);
template void StateFields::set_tensor(char const* name, apf::MeshEntity* e, unsigned n, Intrepid2::Tensor<FadType> const& v);

template void StateFields::get_scalar(char const* name, apf::MeshEntity* e, unsigned n, double& v);
template void StateFields::get_scalar(char const* name, apf::MeshEntity* e, unsigned n, FadType& v);
template void StateFields::get_vector(char const* name, apf::MeshEntity* e, unsigned n, Intrepid2::Vector<double>& v);
template void StateFields::get_vector(char const* name, apf::MeshEntity* e, unsigned n, Intrepid2::Vector<FadType>& v);
template void StateFields::get_tensor(char const* name, apf::MeshEntity* e, unsigned n, Intrepid2::Tensor<double>& v);
template void StateFields::get_tensor(char const* name, apf::MeshEntity* e, unsigned n, Intrepid2::Tensor<FadType>& v);

}
