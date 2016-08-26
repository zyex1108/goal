#include "size_field.hpp"
#include "mesh.hpp"
#include "assert_param.hpp"
#include "control.hpp"

#include <ma.h>
#include <spr.h>

namespace goal {

static RCP<ParameterList> get_valid_uniform_size_params()
{
  RCP<ParameterList> p = rcp(new ParameterList);
  p->set<std::string>("type", "");
  return p;
}

static void validate_uniform_size_params(
    RCP<const ParameterList> p)
{
  p->validateParameters(*get_valid_uniform_size_params(), 0);
}

static RCP<ParameterList> get_valid_spr_size_params()
{
  RCP<ParameterList> p = rcp(new ParameterList);
  p->set<std::string>("type", "");
  p->set<std::string>("variable", "");
  p->set<unsigned>("target elements", 0);
  return p;
}

static void validate_spr_size_params(
    RCP<const ParameterList> p,
    RCP<Mesh> m)
{
  assert_param(p, "variable");
  assert_param(p, "target elements");
  std::string const& var = p->get<std::string>("variable");
  CHECK(m->get_apf_mesh()->findField(var.c_str()));
  p->validateParameters(*get_valid_spr_size_params(), 0);
}

void validate_size_params(
    RCP<const ParameterList> p,
    RCP<Mesh> m)
{
  std::string const& type = p->get<std::string>("type");
  if (type == "uniform")
    validate_uniform_size_params(p);
  else if (type == "spr")
    validate_spr_size_params(p, m);
  else
    fail("unkown size field: %s", type.c_str());
}

apf::Field* get_spr_size_field(
    RCP<const ParameterList> p,
    RCP<Mesh> m)
{
  std::string const& var = p->get<std::string>("variable");
  apf::Field* f = m->get_apf_mesh()->findField(var.c_str());
  CHECK(f);
  unsigned t = p->get<unsigned>("target elements");
  return spr::getTargetSPRSizeField(f, t);
}

apf::Field* get_size_field(
    RCP<const ParameterList> p,
    RCP<Mesh> m)
{
  RCP<const ParameterList> sp = rcpFromRef(p->sublist("size field"));
  std::string const& type = sp->get<std::string>("type");
  if (type == "spr")
    return get_spr_size_field(sp, m);
  else
    fail("unkown size field: %s", type.c_str());
  return 0;
}

}
