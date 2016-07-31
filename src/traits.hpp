#ifndef goal_traits_hpp
#define goal_traits_hpp

#include "dimension.hpp"

#include <Sacado.hpp>
#include <Sacado_mpl_vector.hpp>
#include <Sacado_mpl_find.hpp>
#include <Phalanx_config.hpp>
#include <Phalanx_Traits.hpp>

namespace goal {

struct Workset;

struct GoalTraits : public PHX::TraitsBase
{
  typedef double RealType;
  typedef Sacado::Fad::SLFad<RealType, GOAL_FAD_SIZE> FadType;
  struct Forward {typedef RealType ScalarT;};
  struct Derivative {typedef FadType ScalarT;};
  typedef Sacado::mpl::vector<Forward, Derivative> EvalTypes;
  typedef void* SetupData;
  typedef Workset& PreEvalData;
  typedef Workset& PostEvalData;
  typedef Workset& EvalData;
};

}

namespace PHX {

template <>
struct eval_scalar_types<goal::GoalTraits::Forward>
{
  typedef Sacado::mpl::vector<
    goal::GoalTraits::RealType> type;
};

template <>
struct eval_scalar_types<goal::GoalTraits::Derivative>
{
  typedef Sacado::mpl::vector<
    goal::GoalTraits::FadType,
    goal::GoalTraits::RealType> type;
};

}

#define GOAL_INSTANTIATE_FORWARD(name) \
  template class name<goal::GoalTraits::Forward, goal::GoalTraits>;

#define GOAL_INSTANTIATE_DERIVATIVE(name) \
  template class name<goal::GoalTraits::Derivative, goal::GoalTraits>;

#define GOAL_INSTANTIATE_ALL(name) \
  GOAL_INSTANTIATE_FORWARD(name) \
  GOAL_INSTANTIATE_DERIVATIVE(name)

#endif
