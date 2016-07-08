#include "control.hpp"
#include "expression.hpp"

#include <PCU.h>
#include <cstdlib>
#include <cstdarg>

namespace goal {

void initialize()
{
  MPI_Init(0,0);
  PCU_Comm_Init();
  expression_init();
}

void finalize()
{
  PCU_Comm_Free();
  MPI_Finalize();
}

void print(const char* format, ...)
{
  if (PCU_Comm_Self())
    return void();
  va_list ap;
  va_start(ap, format);
  vfprintf(stdout, format, ap);
  va_end(ap);
  printf("\n");
}

void fail(const char* format, ...)
{
  va_list ap;
  va_start(ap, format);
  vfprintf(stderr, format, ap);
  va_end(ap);
  printf("\n");
  abort();
}

double time()
{
  return PCU_Time();
}

}
