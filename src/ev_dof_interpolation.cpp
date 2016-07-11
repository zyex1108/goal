#include "ev_dof_interpolation.hpp"
#include "traits.hpp"
#include "workset.hpp"
#include "layouts.hpp"
#include "phx_utils.hpp"

namespace goal {

const std::string dof_names[3] =
{ "DOFs",
  "DOFs dot",
  "DOFs dot dot"};

PHX_EVALUATOR_CTOR(DOFInterpolation, p) :
  dl      (p.get<RCP<Layouts> >("Layouts")),
  names   (p.get<Teuchos::Array<std::string> >("DOF Names")),
  index   (p.get<unsigned>("Sol Index")),
  BF      (p.get<std::string>("BF Name"), dl->node_qp_scalar)
{
  num_nodes = dl->node_qp_vector->dimension(1);
  num_qps = dl->node_qp_vector->dimension(2);
  num_dims = dl->node_qp_vector->dimension(3);
  num_eqs = names.size();

  get_grad_field<double>(BF, dl, gBF);
  this->addDependentField(BF);
  this->addDependentField(gBF);

  nodal.resize(num_eqs);
  dofs.resize(num_eqs);
  gdofs.resize(num_eqs);

  for (unsigned i=0; i < num_eqs; ++i) {
    get_field<ScalarT>(names[i], dl, nodal[i]);
    get_field<ScalarT>(names[i], dl, dofs[i]);
    get_grad_field<ScalarT>(names[i], dl, gdofs[i]);
    this->addDependentField(nodal[i]);
    this->addEvaluatedField(dofs[i]);
    this->addEvaluatedField(gdofs[i]);
  }

  this->setName("Interpolate " + dof_names[index]);
}

PHX_POST_REGISTRATION_SETUP(DOFInterpolation, data, fm)
{
  this->utils.setFieldData(BF, fm);
  this->utils.setFieldData(gBF, fm);
  for (unsigned i=0; i < num_eqs; ++i) {
    this->utils.setFieldData(nodal[i], fm);
    this->utils.setFieldData(dofs[i], fm);
    this->utils.setFieldData(gdofs[i], fm);
  }
}

PHX_EVALUATE_FIELDS(DOFInterpolation, workset)
{
  for (unsigned elem=0; elem < workset.size; ++elem) {
  for (unsigned qp=0; qp < num_qps; ++qp) {
  for (unsigned eq=0; eq < num_eqs; ++eq) {
    dofs[eq](elem,qp) = nodal[eq](elem,0) * BF(elem,0,qp);
    for (unsigned node=1; node < num_nodes; ++node)
      dofs[eq](elem, qp) += nodal[eq](elem, node) * BF(elem,node,qp);
    for (unsigned dim=0; dim < num_dims; ++dim) {
      gdofs[eq](elem,qp,dim) = nodal[eq](elem,0) * gBF(elem,0,qp,dim);
      for (unsigned node=1; node < num_nodes; ++node)
        gdofs[eq](elem,qp,dim) += nodal[eq](elem,node) * gBF(elem,node,qp,dim);
  }}}}
}

GOAL_INSTANTIATE_ALL(DOFInterpolation)

}
