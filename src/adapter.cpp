#include "adapter.hpp"
#include "mesh.hpp"
#include "mechanics.hpp"
#include "solution_info.hpp"
#include "solution_attachment.hpp"
#include "size_field.hpp"
#include "workset.hpp"
#include "assert_param.hpp"
#include "control.hpp"

#include <ma.h>

namespace goal {

static RCP<ParameterList> get_valid_params()
{
  RCP<ParameterList> p = rcp(new ParameterList);
  Teuchos::Array<std::string> dummy(0);
  p->sublist("size field");
  p->set<unsigned>("max iters", 0);
  p->set<Teuchos::Array<std::string> >("lb", dummy);
  return p;
}

static void validate_params(
    RCP<const ParameterList> p,
    RCP<Mesh> mesh)
{
  assert_sublist(p, "size field");
  assert_param(p, "lb");
  Teuchos::Array<std::string> lb;
  lb = p->get<Teuchos::Array<std::string> >("lb");
  CHECK(lb.size() == 3);
  for (unsigned i=0; i < 3; ++i)
    CHECK((lb[i] == "none") || (lb[i] == "parma") || (lb[i] == "zoltan"));
  RCP<const ParameterList> sp = rcpFromRef(p->sublist("size field"));
  validate_size_params(sp, mesh);
  p->validateParameters(*get_valid_params(), 0);
}

Adapter::Adapter(
    RCP<const ParameterList> p,
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s) :
  params(p),
  mesh(m),
  mechanics(mech),
  sol_info(s),
  max_iters(0),
  size_field(0)
{
  validate_params(params, mesh);
  max_iters = params->get<unsigned>("max iters");
  size_field_type = params->sublist("size field").get<std::string>("type");
  load_balance = params->get<Teuchos::Array<std::string> >("lb");
}

void Adapter::pre_adapt()
{
  AttachInfo info = {mesh, mechanics, sol_info};
  attach_solutions_to_shape(info);
}

ma::Input* Adapter::create_input()
{
  ma::Input* in;
  if (size_field_type == "uniform")
    in = ma::configureUniformRefine(mesh->get_apf_mesh());
  else {
    size_field = get_size_field(params, mesh);
    in = ma::configure(mesh->get_apf_mesh(), size_field);
  }
  in->shouldRunPreZoltan = ("zoltan" == load_balance[0]);
  in->shouldRunPreParma = ("parma" == load_balance[0]);
  in->shouldRunMidZoltan = ("zoltan" == load_balance[1]);
  in->shouldRunMidParma = ("parma" == load_balance[1]);
  in->shouldRunPostZoltan = ("zoltan" == load_balance[2]);
  in->shouldRunPostParma = ("parma" == load_balance[2]);
  return in;
}

void Adapter::adapt(unsigned step_number)
{
  pre_adapt();
  ma::Input* in = create_input();
  ma::adapt(in);
  post_adapt();
}

void Adapter::post_adapt()
{
  apf::destroyField(size_field);
  mesh->update();
  sol_info->resize(mesh, false);
  AttachInfo info = {mesh, mechanics, sol_info};
  fill_solutions_from_fields(info);
  remove_solutions_from_mesh(info);
}

RCP<Adapter> adapter_create(
    RCP<const ParameterList> p,
    RCP<Mesh> m,
    RCP<Mechanics> mech,
    RCP<SolutionInfo> s)
{
  RCP<const ParameterList> ap = rcpFromRef(p->sublist("adapt"));
  return rcp(new Adapter(ap, m, mech, s));
}

}
