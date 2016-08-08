#ifndef goal_solution_attachment_hpp
#define goal_solution_attachment_hpp

#include <Teuchos_RCP.hpp>

namespace goal {

using Teuchos::RCP;

class Mesh;
class Mechanics;
class SolutionInfo;

struct AttachInfo
{
  RCP<Mesh> mesh;
  RCP<Mechanics> mech;
  RCP<SolutionInfo> sol_info;
};

void attach_solutions_to_mesh(AttachInfo& i);
void attach_dual_solutions_to_mesh(AttachInfo& i);

void attach_solutions_to_shape(AttachInfo& i);
void attach_dual_solutions_to_shape(AttachInfo& i);

void fill_solution_vectors(AttachInfo& i);
void fill_dual_solution_vectors(AttachInfo& i);

void remove_solutions_from_mesh(AttachInfo& i);
void remove_dual_solutions_from_mesh(AttachInfo& i);

}

#endif
