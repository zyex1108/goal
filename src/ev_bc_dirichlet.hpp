#ifndef goal_ev_bc_dirichlet_hpp
#define goal_ev_bc_dirichlet_hpp

#include "phx_macros.hpp"
#include "traits.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class Mesh;
class Mechanics;
struct Layouts;
template <typename EvalT, typename Traits> class BCDirichlet;

template <typename Traits>
class BCDirichlet<GoalTraits::Forward, Traits> :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<GoalTraits::Forward, Traits>
{
  public:

    BCDirichlet(ParameterList const& p);

    void postRegistrationSetup(
        typename Traits::SetupData d,
        PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);

  private:

    typedef typename GoalTraits::Forward::ScalarT ScalarT;

    RCP<Layouts> dl;
    RCP<Mesh> mesh;
    RCP<Mechanics> mechanics;
    RCP<const ParameterList> params;

    void validate_params();

    void apply_bc(
        typename Traits::EvalData d,
        Teuchos::Array<std::string> const& a);
};

template <typename Traits>
class BCDirichlet<GoalTraits::Derivative, Traits> :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<GoalTraits::Derivative, Traits>
{
  public:

    BCDirichlet(ParameterList const& p);

    void postRegistrationSetup(
        typename Traits::SetupData d,
        PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);

  private:

    typedef typename GoalTraits::Derivative::ScalarT ScalarT;

    RCP<Layouts> dl;
    RCP<Mesh> mesh;
    RCP<Mechanics> mechanics;
    RCP<const ParameterList> params;

    void validate_params();

    void apply_bc(
        typename Traits::EvalData d,
        Teuchos::Array<std::string> const& a);

};

}

#endif
