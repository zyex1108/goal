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
  names   (p.get<Teuchos::Array<std::string> >("Dual Names")),
  BF      (p.get<std::string>("BF Name"), dl->node_qp_scalar)
{
  num_nodes = dl->node_qp_scalar->dimension(1);
  num_qps = dl->node_qp_vector->dimension(2);
  num_dims = dl->node_qp_vector->dimension(3);
  num_eqs = names.size();

  get_grad_field<double>(BF, dl, gBF);
  this->addDependentField(BF);
  this->addDependentField(gBF);

  nodal.resize(num_eqs);
  duals.resize(num_eqs);
  gduals.resize(num_eqs);

  for (unsigned i=0; i < num_eqs; ++i) {
    get_field<ScalarT>(names[i], dl, nodal[i]);
    get_field<ScalarT>(names[i], dl, duals[i]);
    get_grad_field<ScalarT>(names[i], dl, gduals[i]);
    this->addEvaluatedField(nodal[i]);
    this->addEvaluatedField(duals[i]);
    this->addEvaluatedField(gduals[i]);
  }

  this->setName("Gather Dual");
}

PHX_POST_REGISTRATION_SETUP(GatherDual, data, fm)
{
  this->utils.setFieldData(BF, fm);
  this->utils.setFieldData(gBF, fm);
  for (unsigned i=0; i < num_eqs; ++i) {
    this->utils.setFieldData(nodal[i], fm);
    this->utils.setFieldData(duals[i], fm);
    this->utils.setFieldData(gduals[i], fm);
  }
}

PHX_EVALUATE_FIELDS(GatherDual, workset)
{
  CHECK(workset.z != Teuchos::null);
  ArrayRCP<const ST> dual = workset.z->get1dView();
  CHECK(dual != Teuchos::null);

  for (unsigned elem=0; elem < workset.size; ++elem) {
    apf::MeshEntity* e = workset.ents[elem];
    for (unsigned node=0; node < num_nodes; ++node) {
      for (unsigned eq=0; eq < num_eqs; ++eq) {
        LO lid = mesh->get_lid(e, node, eq);
        nodal[eq](elem, node) = dual[lid];
      }
    }
  }

  for (unsigned elem=0; elem < workset.size; ++elem) {
  for (unsigned qp=0; qp < num_qps; ++qp) {
  for (unsigned eq=0; eq < num_eqs; ++eq) {
    duals[eq](elem,qp) = nodal[eq](elem,0)*BF(elem,0,qp);
    for (unsigned node=1; node < num_nodes; ++node)
      duals[eq](elem,qp) += nodal[eq](elem,node)*BF(elem,node,qp);
    for (unsigned dim=0; dim < num_dims; ++dim) {
      gduals[eq](elem,qp,dim) = nodal[eq](elem,0)*gBF(elem,0,qp,dim);
      for (unsigned node=1; node < num_nodes; ++node)
        gduals[eq](elem,qp,dim) += nodal[eq](elem,node)*gBF(elem,node,qp,dim);
  }}}}
}

GOAL_INSTANTIATE_ALL(GatherDual)

}
