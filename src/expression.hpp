#ifndef goal_expression_hpp
#define goal_expression_hpp

#include <string>

namespace goal {

void expression_init();

double expression_eval(
    std::string const& val,
    const double x,
    const double y,
    const double z,
    const double t);

}

#endif
