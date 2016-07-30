#ifndef goal_mechanics_hpp
#define goal_mechanics_hpp

#include "traits.hpp"
#include <Phalanx_FieldManager.hpp>

namespace goal {

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::ParameterList;

class Mesh;
class StateFields;

typedef RCP<PHX::FieldManager<GoalTraits> > FieldManager;
typedef ArrayRCP<FieldManager> FieldManagers;

class Mechanics
{
  public:

    Mechanics(
        RCP<const ParameterList> p,
        RCP<Mesh> m,
        bool supports_dynamics);

    unsigned get_num_eqs();
    
    void build_primal();
    void build_dual();
    void build_error();

    void project_state();
    void update_state();

    Teuchos::Array<std::string> const& get_dof_names();
    Teuchos::Array<std::string> const& get_var_names(unsigned sol_idx);
    unsigned get_offset(std::string const& var_name);

    FieldManagers get_volumetric() {return vfms;}
    FieldManager get_dirichlet() {return dfm;}
    FieldManager get_neumann() {return nfm;}

  private:

    RCP<const ParameterList> params;
    RCP<Mesh> mesh;

    bool supports_dynamics;
    bool have_pressure_eq;
    bool have_temperature;
    bool small_strain;

    bool is_primal;
    bool is_dual;
    bool is_error;

    unsigned num_eqs;

    Teuchos::Array<std::string> var_names[3];
    std::map<std::string, unsigned> offsets;
    std::map<std::string, Teuchos::Array<std::string> > fields;

    std::string model;
    Teuchos::RCP<StateFields> state_fields;

    FieldManagers vfms;
    FieldManager nfm;
    FieldManager dfm;

    void set_primal();
    void set_dual();
    void set_error();

    void setup_params();
    void setup_variables();
    void setup_fields();
    void setup_states();

    template <typename EvalT>
    void register_volumetric(std::string const& set, FieldManager fm);

    template <typename EvalT>
    void register_neumann(FieldManager fm);

    template <typename EvalT>
    void register_dirichlet(FieldManager fm);

};

RCP<Mechanics> mechanics_create(
    RCP<const ParameterList> p,
    RCP<Mesh> m,
    bool supports_dynamics);

}

#endif
