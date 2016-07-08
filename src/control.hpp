#ifndef goal_control_hpp
#define goal_control_hpp

#include "config.hpp"

namespace goal {

void initialize();
void finalize();
void print(char const* format, ...);
void fail(char const* format, ...) __attribute__((noreturn));
double time();

}

#define ALWAYS_CHECK(c) \
  ((c) ? ((void)0) :    \
   goal::fail("assertion (%s) failed at %s:%d\n",#c,__FILE__,__LINE__))

#ifdef GOAL_DISABLE_CHECKS
#define CHECK(c) ((void)0)
#else
#define CHECK(c) (ALWAYS_CHECK(c))
#endif

#endif
