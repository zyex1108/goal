#include "mechanics.hpp"
#include "mesh.hpp"
#include "state_fields.hpp"
#include "control.hpp"
#include "assert_param.hpp"

namespace goal {

static RCP<ParameterList> get_valid_params(RCP<Mesh> m)
{
  RCP<ParameterList> p = rcp(new ParameterList);
  p->set<std::string>("model", "");
  p->sublist("dirichlet bcs");
  p->sublist("neumann bcs");
  for (unsigned i=0; i < m->get_num_elem_sets(); ++i) {
    std::string const& set  = m->get_elem_set_name(i);
    p->sublist(set);
  }
  return p;
}

static void validate_params(RCP<const ParameterList> p, RCP<Mesh> m)
{
  assert_param(p, "model");
  assert_param(p, "dirichlet bcs");
  for (unsigned i=0; i < m->get_num_elem_sets(); ++i)
    assert_param(p, m->get_elem_set_name(i));
  p->validateParameters(*get_valid_params(m), 0);
}

Mechanics::Mechanics(
    RCP<const ParameterList> p,
    RCP<Mesh> m,
    bool support) :
  params(p),
  mesh(m),
  supports_dynamics(support)
{
  validate_params(params, mesh);
  setup_dofs();
  setup_fields();
  setup_states();
  mesh->set_num_eqs(num_eqs);
  mesh->update();
}

unsigned Mechanics::get_num_eqs()
{
  return dof_names.size();
}

void Mechanics::set_primal()
{
  is_primal = true;
  is_dual = false;
  is_error = false;
}

void Mechanics::set_dual()
{
  is_primal = false;
  is_dual = false;
  is_error = false;
}

void Mechanics::set_error()
{
  is_primal = false;
  is_dual = false;
  is_error = true;
}

void Mechanics::build_primal()
{
  print("building primal pde fields");
}

void Mechanics::build_dual()
{
  print("building dual pde fields");
}

void Mechanics::build_error()
{
  print("building error fields");
}

void Mechanics::project_state()
{
}

void Mechanics::update_state()
{
}

void Mechanics::setup_dofs()
{
  unsigned d = mesh->get_num_dims();
  dof_names.push_back("ux");
  if (d > 1) dof_names.push_back("uy");
  if (d > 2) dof_names.push_back("uz");
  dof_names.push_back("p");
  for (unsigned i=0; i < dof_names.size(); ++i)
    dof_offsets[dof_names[i]] = i;
  num_eqs = dof_names.size();
}

void Mechanics::setup_fields()
{
  Teuchos::Array<std::string> disp;
  for (unsigned i=0; i < mesh->get_num_dims(); ++i)
    disp.push_back(dof_names[i]);
  fields["disp"] = disp;
  Teuchos::Array<std::string> pressure;
  pressure.push_back("p");
  fields["pressure"] = pressure;
}

void Mechanics::setup_states()
{
  model = params->get<std::string>("model");
  state_fields = rcp(new StateFields(mesh));
  if (model == "linear elastic")
    state_fields->add("cauchy", TENSOR, true);
  else if (model == "J2") {
    state_fields->add("cauchy", TENSOR, true);
    state_fields->add("Fp", TENSOR, true, true);
    state_fields->add("eqps", SCALAR, true);
  }
}

RCP<Mechanics> mechanics_create(
    RCP<const ParameterList> p,
    RCP<Mesh> m,
    bool supports_dynamics)
{
  RCP<const ParameterList> mp = rcpFromRef(p->sublist("mechanics"));
  return rcp(new Mechanics(mp, m, supports_dynamics));
}

}
