#include "adapter.hpp"
#include "mesh.hpp"
#include "mechanics.hpp"
#include "solution_info.hpp"
#include "solution_attachment.hpp"
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

static void validate_params(RCP<const ParameterList> p)
{
  assert_sublist(p, "size field");
  assert_param(p, "lb");
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
  max_iters(0)
{
  validate_params(params);
  max_iters = params->get<unsigned>("max iters");
}

void Adapter::pre_adapt()
{
  AttachInfo info = {mesh, mechanics, sol_info};
  attach_solutions_to_shape(info);
}

void Adapter::adapt(unsigned step_number)
{
  pre_adapt();
  ma::runUniformRefinement(mesh->get_apf_mesh());
  post_adapt();
}

void Adapter::post_adapt()
{
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
