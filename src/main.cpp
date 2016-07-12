#include "control.hpp"
#include "solver.hpp"
#include <Teuchos_XMLParameterListHelpers.hpp>

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::ParameterList;

int main(int argc, char** argv)
{
  goal::initialize();
  goal::print("Goal!");
  ALWAYS_CHECK(argc == 2);
  int status = EXIT_SUCCESS;
  try {
    char const* in = argv[1];
    RCP<ParameterList> p = rcp(new ParameterList);
    goal::print("reading input file: %s", in);
    Teuchos::updateParametersFromXmlFile(in, p.ptr());
    RCP<goal::Solver> solver = goal::solver_create(p);
    solver->solve();
  }
  catch (std::exception const& ex) {
    goal::print("caught exception:");
    goal::print("%s", ex.what());
    status = EXIT_FAILURE;
  }
  catch (...) {
    goal::print("caught unknown exception");
    status = EXIT_FAILURE;
  }
  goal::finalize();
  return status;
}
