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
  p->sublist("qoi");
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
  small_strain(false)
{
  validate_params(params, mesh);
  setup_params();
  setup_variables();
  setup_fields();
  setup_states();
  mesh->set_num_eqs(num_eqs);
  mesh->update();
}

unsigned Mechanics::get_num_eqs()
{
  return num_eqs;
}

void Mechanics::build_primal()
{
  double t0 = time();
  set_primal();
  typedef GoalTraits::Forward F;
  typedef GoalTraits::Derivative D;
  vfms.resize(mesh->get_num_elem_sets());
  for (unsigned i=0; i < mesh->get_num_elem_sets(); ++i) {
    vfms[i] = rcp(new PHX::FieldManager<GoalTraits>);
    std::string const& set = mesh->get_elem_set_name(i);
    register_volumetric<F>(set, vfms[i]);
    register_volumetric<D>(set, vfms[i]);
  }
  nfm = rcp(new PHX::FieldManager<GoalTraits>);
  dfm = rcp(new PHX::FieldManager<GoalTraits>);
  register_neumann<F>(nfm);
  register_neumann<D>(nfm);
  register_dirichlet<F>(dfm);
  register_dirichlet<D>(dfm);
  double t1 = time();
  print("primal pde fields built in %f seconds", t1-t0);
}

void Mechanics::build_dual()
{
  double t0 = time();
  set_dual();
  typedef GoalTraits::Derivative D;
  vfms.resize(mesh->get_num_elem_sets());
  for (unsigned i=0; i < mesh->get_num_elem_sets(); ++i) {
    vfms[i] = rcp(new PHX::FieldManager<GoalTraits>);
    std::string const& set = mesh->get_elem_set_name(i);
    register_volumetric<D>(set, vfms[i]);
  }
  dfm = rcp(new PHX::FieldManager<GoalTraits>);
  register_dirichlet<D>(dfm);
  double t1 = time();
  print("dual pde fields built in %f seconds", t1-t0);
}

void Mechanics::build_error()
{
  double t0 = time();
  double t1 = time();
  print("error fields built in %f seconds", t1-t0);
}

void Mechanics::project_state()
{
}

void Mechanics::update_state()
{
  state_fields->update();
}

Teuchos::Array<std::string> const& Mechanics::get_dof_names()
{
  return var_names[0];
}

Teuchos::Array<std::string> const& Mechanics::get_var_names(unsigned i)
{
  if (i > 0) CHECK(supports_dynamics);
  CHECK((i>=0) && (i<=2));
  return var_names[i];
}

unsigned Mechanics::get_offset(std::string const& var_name)
{
  CHECK(offsets.count(var_name));
  return offsets[var_name];
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
  is_dual = true;
  is_error = false;
}

void Mechanics::set_error()
{
  is_primal = false;
  is_dual = false;
  is_error = true;
}

void Mechanics::setup_params()
{
  model = params->get<std::string>("model");
  if (params->isParameter("mixed formulation"))
    have_pressure_eq = params->get<bool>("mixed formulation");
  if (params->isParameter("temperature"))
    have_temperature = true;
}

void Mechanics::setup_variables()
{
  unsigned d = mesh->get_num_dims();

  var_names[0].push_back("ux");
  if (d > 1) var_names[0].push_back("uy");
  if (d > 2) var_names[0].push_back("uz");
  if (supports_dynamics) {
    var_names[1].push_back("vx");
    if (d > 1) var_names[1].push_back("vy");
    if (d > 2) var_names[1].push_back("vz");
    var_names[2].push_back("ax");
    if (d > 1) var_names[2].push_back("ay");
    if (d > 2) var_names[2].push_back("az");
  }

  if (have_pressure_eq) {
    var_names[0].push_back("p");
    if (supports_dynamics) {
      var_names[1].push_back("dp");
      var_names[2].push_back("ddp");
    }
  }

  if (supports_dynamics) {
    CHECK(var_names[0].size() == var_names[1].size());
    CHECK(var_names[1].size() == var_names[2].size());
  }

  for (unsigned i=0; i < var_names[0].size(); ++i)
    offsets[var_names[0][i]] = i;
  for (unsigned i=0; i < var_names[1].size(); ++i)
    offsets[var_names[1][i]] = i;
  for (unsigned i=0; i < var_names[2].size(); ++i)
    offsets[var_names[2][i]] = i;

  num_eqs = var_names[0].size();
}

void Mechanics::setup_fields()
{
  unsigned d = mesh->get_num_dims();

  Teuchos::Array<std::string> disp;
  disp.push_back("ux");
  if (d > 1) disp.push_back("uy");
  if (d > 2) disp.push_back("uz");
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
    state_fields->add("cauchy", TENSOR, false);
  else if (model == "j2") {
    state_fields->add("cauchy", TENSOR, false);
    state_fields->add("Fp", TENSOR, true, true);
    state_fields->add("eqps", SCALAR, true);
  }
  else if (model == "creep") {
    state_fields->add("cauchy", TENSOR, false);
    state_fields->add("Fp", TENSOR, true, true);
    state_fields->add("eqps", SCALAR, true);
  }
  else
    fail("unknown model: %s", model.c_str());
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
