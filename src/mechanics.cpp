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
  p->set<bool>("mixed formulation", false);
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
  supports_dynamics(support),
  have_pressure_eq(false),
  have_temperature(false),
  small_strain(true)
{
  validate_params(params, mesh);
  setup_params();
  setup_dofs();
  setup_fields();
  setup_states();
  mesh->set_num_eqs(num_eqs);
  mesh->update();
}

unsigned Mechanics::get_num_eqs()
{
  return num_eqs;
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
  set_primal();
  typedef GoalTraits::Residual R;
  typedef GoalTraits::Jacobian J;
  pfms.resize(mesh->get_num_elem_sets());
  for (unsigned i=0; i < mesh->get_num_elem_sets(); ++i) {
    pfms[i] = rcp(new PHX::FieldManager<GoalTraits>);
    std::string const& set = mesh->get_elem_set_name(i);
    register_volumetric<R>(set, pfms[i]);
    register_volumetric<J>(set, pfms[i]);
  }
  nfm = rcp(new PHX::FieldManager<GoalTraits>);
  dfm = rcp(new PHX::FieldManager<GoalTraits>);
  register_neumann<R>(nfm);
  register_neumann<J>(nfm);
  register_dirichlet<R>(dfm);
  register_dirichlet<J>(dfm);
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

Teuchos::Array<std::string> const& Mechanics::get_dof_names()
{
  return dof_names;
}

Teuchos::Array<std::string> const& Mechanics::get_dof_dot_names()
{
  return dof_dot_names;
}

Teuchos::Array<std::string> const& Mechanics::get_dof_dot_dot_names()
{
  return dof_dot_dot_names;
}

void Mechanics::setup_params()
{
  model = params->get<std::string>("model");
  if (params->isParameter("mixed formulation"))
    have_pressure_eq = params->get<bool>("mixed formulation");
  if (params->isParameter("temperature"))
    have_temperature = true;
}

void Mechanics::setup_dofs()
{
  unsigned d = mesh->get_num_dims();

  dof_names.push_back("ux");
  if (d > 1) dof_names.push_back("uy");
  if (d > 2) dof_names.push_back("uz");
  if (supports_dynamics) {
    dof_dot_names.push_back("vx");
    if (d > 1) dof_dot_names.push_back("vy");
    if (d > 2) dof_dot_names.push_back("vz");
    dof_dot_dot_names.push_back("ax");
    if (d > 1) dof_dot_dot_names.push_back("ay");
    if (d > 2) dof_dot_dot_names.push_back("az");
  }

  if (have_pressure_eq) dof_names.push_back("p");
  if (supports_dynamics) {
    dof_dot_names.push_back("p_dot");
    dof_dot_dot_names.push_back("p_dot_dot");
  }

  for (unsigned i=0; i < dof_names.size(); ++i)
    dof_offsets[dof_names[i]] = i;

  num_eqs = dof_names.size();
}

void Mechanics::setup_fields()
{
  unsigned d = mesh->get_num_dims();

  Teuchos::Array<std::string> disp;
  disp.push_back("ux");
  if (d > 1) disp.push_back("uy");
  if (d > 1) disp.push_back("uz");
  fields["disp"] = disp;

  if (supports_dynamics) {
    Teuchos::Array<std::string> vel;
    vel.push_back("vx");
    if (d > 1) vel.push_back("vy");
    if (d > 2) vel.push_back("vz");
    fields["vel"] = vel;

    Teuchos::Array<std::string> acc;
    acc.push_back("ax");
    if (d > 1) acc.push_back("ay");
    if (d > 2) acc.push_back("az");
    fields["acc"] = acc;
  }

  if (have_pressure_eq) {
    Teuchos::Array<std::string> p;
    p.push_back("p");
    fields["p"] = p;
  }
}

void Mechanics::setup_states()
{
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
