#ifndef goal_ev_bc_traction_hpp
#define goal_ev_bc_traction_hpp

#include "phx_macros.hpp"

#include <apf.h>

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

struct Layouts;
class Mesh;
class Mechanics;

PHX_EVALUATOR_CLASS(BCTraction)

  private:

    RCP<Layouts> dl;
    RCP<Mesh> mesh;
    RCP<Mechanics> mechanics;
    RCP<const ParameterList> params;

    unsigned num_dims;

    apf::Mesh* apf_mesh;
    apf::Vector3 local;
    apf::Vector3 global;
    apf::Vector3 traction;
    apf::NewArray<long> numbers;
    apf::NewArray<double> BF;

    void validate_params();
    void apply_bc(
        typename Traits::EvalData d,
        Teuchos::Array<std::string> const& a);

PHX_EVALUATOR_CLASS_END

}

#endif
