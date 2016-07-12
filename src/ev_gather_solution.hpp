#ifndef goal_ev_gather_solution_hpp
#define goal_ev_gather_solution_hpp

#include "phx_macros.hpp"
#include "traits.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::ParameterList;

class Mesh;
struct Layouts;
template <typename EvalT, typename Traits> class GatherSolution;

/* residual specialization */

template <typename Traits>
class GatherSolution<GoalTraits::Residual, Traits> :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<GoalTraits::Residual, Traits>
{
  public:

    GatherSolution(ParameterList const& p);

    void postRegistrationSetup(
        typename Traits::SetupData d,
        PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);

  private:

    typedef typename GoalTraits::Residual::ScalarT ScalarT;

    RCP<Layouts> dl;
    RCP<Mesh> mesh;
    Teuchos::Array<std::string> names;
    unsigned index;

    unsigned num_nodes;
    unsigned num_eqs;

    std::vector<PHX::MDField<ScalarT, Elem, Node> > u;
};


/* jacobian specialization */

template <typename Traits>
class GatherSolution<GoalTraits::Jacobian, Traits> :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<GoalTraits::Jacobian, Traits>
{
  public:

    GatherSolution(ParameterList const& p);

    void postRegistrationSetup(
        typename Traits::SetupData d,
        PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);

  private:

    typedef typename GoalTraits::Jacobian::ScalarT ScalarT;

    RCP<Layouts> dl;
    RCP<Mesh> mesh;
    Teuchos::Array<std::string> names;

    unsigned num_nodes;
    unsigned num_eqs;
    unsigned index;

    std::vector<PHX::MDField<ScalarT, Elem, Node> > u;
};

}

#endif
