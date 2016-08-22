#include "adapter.hpp"
#include "mesh.hpp"
#include "mechanics.hpp"
#include "solution_info.hpp"
#include "workset.hpp"
#include "assert_param.hpp"
#include "control.hpp"

namespace goal {

static RCP<ParameterList> get_valid_params()
{
  RCP<ParameterList> p = rcp(new ParameterList);
  Teuchos::Array<std::string> dummy(0);
  p->sublist("size field");
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
  sol_info(s)
{
  validate_params(params);
}

void Adapter::adapt(unsigned step_number)
{
  std::cout << "I am adapting" << std::endl;
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
