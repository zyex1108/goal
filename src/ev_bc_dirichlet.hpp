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
class BCDirichlet<GoalTraits::Residual, Traits> :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<GoalTraits::Residual, Traits>
{
  public:

    BCDirichlet(ParameterList const& p);

    void postRegistrationSetup(
        typename Traits::SetupData d,
        PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);

  private:

    typedef typename GoalTraits::Residual::ScalarT ScalarT;

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
class BCDirichlet<GoalTraits::Jacobian, Traits> :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<GoalTraits::Jacobian, Traits>
{
  public:

    BCDirichlet(ParameterList const& p);

    void postRegistrationSetup(
        typename Traits::SetupData d,
        PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);

  private:

    typedef typename GoalTraits::Jacobian::ScalarT ScalarT;

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
