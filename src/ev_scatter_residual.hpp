#ifndef goal_scatter_residual_hpp
#define goal_scatter_residual_hpp

#include "phx_macros.hpp"
#include "traits.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class Mesh;
struct Layouts;
template <typename EvalT, typename Traits> class ScatterResidual;

template <typename Traits>
class ScatterResidual<GoalTraits::Forward, Traits> :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<GoalTraits::Forward, Traits>
{
  public:

    ScatterResidual(ParameterList const& p);

    void postRegistrationSetup(
        typename Traits::SetupData d,
        PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);

  private:

    typedef typename GoalTraits::Forward::ScalarT ScalarT;

    RCP<Layouts> dl;
    RCP<Mesh> mesh;
    Teuchos::Array<std::string> dof_names;

    unsigned num_nodes;
    unsigned num_eqs;

    std::vector<PHX::MDField<ScalarT, Elem, Node> > resid;
};

template <typename Traits>
class ScatterResidual<GoalTraits::Derivative, Traits> :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<GoalTraits::Derivative, Traits>
{
  public:

    ScatterResidual(ParameterList const& p);

    void postRegistrationSetup(
        typename Traits::SetupData d,
        PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);

  private:

    typedef typename GoalTraits::Derivative::ScalarT ScalarT;

    RCP<Layouts> dl;
    RCP<Mesh> mesh;
    Teuchos::Array<std::string> dof_names;

    unsigned num_nodes;
    unsigned num_eqs;
    unsigned num_dofs;

    std::vector<PHX::MDField<ScalarT, Elem, Node> > resid;
};

}

#endif
