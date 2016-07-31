#ifndef goal_ev_scatter_qoi_hpp
#define goal_ev_scatter_qoi_hpp

#include "phx_macros.hpp"
#include "traits.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

struct Layouts;
class Mesh;
template <typename EvalT, typename Traits> class ScatterQoI;

template <typename Traits>
class ScatterQoI<GoalTraits::Forward, Traits> :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<GoalTraits::Forward, Traits>
{
  public:

    ScatterQoI(ParameterList const& p);
    
    void postRegistrationSetup(
        typename Traits::SetupData d,
        PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);

  private:

    typedef typename GoalTraits::Forward::ScalarT ScalarT;

    RCP<Layouts> dl;

    PHX::MDField<ScalarT, Elem> qoi;
};

template <typename Traits>
class ScatterQoI<GoalTraits::Derivative, Traits> :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<GoalTraits::Derivative, Traits>
{
  public:

    ScatterQoI(ParameterList const& p);
    
    void postRegistrationSetup(
        typename Traits::SetupData d,
        PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);

  private:

    typedef typename GoalTraits::Derivative::ScalarT ScalarT;

    RCP<Layouts> dl;
    RCP<Mesh> mesh;

    unsigned num_nodes;

    PHX::MDField<ScalarT, Elem> qoi;
};

}

#endif
