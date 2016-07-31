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
class ScatterQoI<GoalTraits::Residual, Traits> :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<GoalTraits::Residual, Traits>
{
  public:

    ScatterQoI(ParameterList const& p);
    
    void postRegistrationSetup(
        typename Traits::SetupData d,
        PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);

  private:

    typedef typename GoalTraits::Residual::ScalarT ScalarT;

    RCP<Layouts> dl;

    PHX::MDField<ScalarT, Elem> qoi;
};

template <typename Traits>
class ScatterQoI<GoalTraits::Jacobian, Traits> :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<GoalTraits::Jacobian, Traits>
{
  public:

    ScatterQoI(ParameterList const& p);
    
    void postRegistrationSetup(
        typename Traits::SetupData d,
        PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);

  private:

    typedef typename GoalTraits::Jacobian::ScalarT ScalarT;

    RCP<Layouts> dl;
    RCP<Mesh> mesh;

    unsigned num_nodes;

    PHX::MDField<ScalarT, Elem> qoi;
};

}

#endif
