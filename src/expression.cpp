#include "expression.hpp"

#include <RTC_FunctionRTC.hh>

namespace goal {

static PG_RuntimeCompiler::Function evaluator(5);

void expression_init()
{
  evaluator.addVar("double", "x");
  evaluator.addVar("double", "y");
  evaluator.addVar("double", "z");
  evaluator.addVar("double", "t");
  evaluator.addVar("double", "val");
}

double expression_eval(
    std::string const& val,
    const double x,
    const double y,
    const double z,
    const double t)
{
  evaluator.addBody(val);
  evaluator.varValueFill(0, x);
  evaluator.varValueFill(1, y);
  evaluator.varValueFill(2, z);
  evaluator.varValueFill(3, t);
  evaluator.varValueFill(4, 0);
  evaluator.execute();
  return evaluator.getValueOfVar("val");
}

}
