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
  struct Residual {typedef RealType ScalarT;};
  struct Jacobian {typedef FadType ScalarT;};
  typedef Sacado::mpl::vector<Residual, Jacobian> EvalTypes;
  typedef void* SetupData;
  typedef Workset& PreEvalData;
  typedef Workset& PostEvalData;
  typedef Workset& EvalData;
};

}

namespace PHX {

template <>
struct eval_scalar_types<goal::GoalTraits::Residual>
{
  typedef Sacado::mpl::vector<
    goal::GoalTraits::RealType> type;
};

template <>
struct eval_scalar_types<goal::GoalTraits::Jacobian>
{
  typedef Sacado::mpl::vector<
    goal::GoalTraits::FadType,
    goal::GoalTraits::RealType> type;
};

}

#define GOAL_INSTANTIATE_RESIDUAL(name) \
  template class name<goal::GoalTraits::Residual, goal::GoalTraits>;

#define GOAL_INSTANTIATE_JACOBIAN(name) \
  template class name<goal::GoalTraits::Jacobian, goal::GoalTraits>;

#define GOAL_INSTANTIATE_ALL(name) \
  GOAL_INSTANTIATE_RESIDUAL(name) \
  GOAL_INSTANTIATE_JACOBIAN(name)

#endif
