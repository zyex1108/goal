#ifndef goal_phx_macros_hpp
#define goal_phx_macros_hpp

#include "dimension.hpp"

#include <Phalanx_config.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>
#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_FieldManager.hpp>

namespace Teuchos {
class ParameterList;
}

/*****************************************************************************/

#define PHX_EVALUATOR_CLASS(NAME) \
  template<typename EvalT, typename Traits> \
  class NAME : \
    public PHX::EvaluatorWithBaseImpl<Traits>, \
    public PHX::EvaluatorDerived<EvalT, Traits> \
  { \
    public: \
      NAME(const Teuchos::ParameterList& p); \
      void postRegistrationSetup( \
          typename Traits::SetupData d, \
          PHX::FieldManager<Traits>& fm); \
      void evaluateFields(typename Traits::EvalData d); \
    private: \
      typedef typename EvalT::ScalarT ScalarT;

/*****************************************************************************/

#define PHX_EVALUATOR_CLASS_END };

/*****************************************************************************/

#define PHX_EVALUATOR_CTOR(NAME,PLIST) \
  template<typename EvalT, typename Traits>  \
  NAME <EvalT, Traits>::NAME(Teuchos::ParameterList const& PLIST)

/*****************************************************************************/

#define PHX_POST_REGISTRATION_SETUP(NAME,SETUP_DATA,FIELD_MANAGER) \
  template<typename EvalT, typename Traits>  \
  void NAME<EvalT, Traits>:: \
  postRegistrationSetup( \
      typename Traits::SetupData SETUP_DATA, \
      PHX::FieldManager<Traits>& FIELD_MANAGER)

/*****************************************************************************/

#define PHX_EVALUATE_FIELDS(NAME,EVAL_DATA) \
  template<typename EvalT, typename Traits> \
  void NAME<EvalT, Traits>:: \
  evaluateFields(typename Traits::EvalData EVAL_DATA)

#endif
