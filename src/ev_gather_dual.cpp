#include "ev_gather_dual.hpp"
#include "mesh.hpp"
#include "workset.hpp"
#include "layouts.hpp"
#include "phx_utils.hpp"
#include "control.hpp"

namespace goal {

PHX_EVALUATOR_CTOR(GatherDual, p) :
  dl      (p.get<RCP<Layouts> >("Layouts")),
  mesh    (p.get<RCP<Mesh> >("Mesh")),
  names   (p.get<Teuchos::Array<std::string> >("Sol Names"))
{
  num_eqs = names.size();
  num_nodes = dl->node_scalar->dimension(1);

  z.resize(num_eqs);
  for (unsigned i=0; i < num_eqs; ++i) {
    get_dual_field(names[i], dl, z[i]);
    this->addEvaluatedField(z[i]);
  }

  this->setName("Gather Dual");
}

PHX_POST_REGISTRATION_SETUP(GatherDual, data, fm)
{
  for (unsigned i=0; i < z.size(); ++i)
    this->utils.setFieldData(z[i], fm);
}

PHX_EVALUATE_FIELDS(GatherDual, workset)
{
  CHECK(workset.z != Teuchos::null);
  ArrayRCP<const ST> dual = z->get1dView();
  CHECK(dual != Teuchos::null);

  for (unsigned elem=0; elem < workset.size; ++elem) {
    apf::MeshEntity* e = workset.ents[elem];
    for (unsigned node=0; node < num_nodes; ++node) {
      for (unsigned eq=0; eq < num_eqs; ++eq) {
        LO lid = mesh->get_lid(e, node, eq);
        z[eq](elem, node) = dual[lid];
      }
    }
  }
}

}
