#include "ev_basis_functions.hpp"
#include "mesh.hpp"
#include "layouts.hpp"
#include "workset.hpp"
#include "traits.hpp"
#include "phx_utils.hpp"

#include <apf.h>
#include <apfMesh2.h>

namespace goal {

PHX_EVALUATOR_CTOR(BasisFunctions, p) :
  dl    (p.get<RCP<Layouts> >("Layouts")),
  mesh  (p.get<RCP<Mesh> >("Mesh")),
  wDv   (p.get<std::string>("Weighted Dv Name"), dl->qp_scalar),
  BF    (p.get<std::string>("BF Name"), dl->node_qp_scalar)
{
  num_nodes = dl->node_qp_vector->dimension(1);
  num_qps = dl->node_qp_vector->dimension(2);
  num_dims = dl->node_qp_vector->dimension(3);

  get_grad_field<double>(BF, dl, gBF);

  p_order = mesh->get_p_order();
  q_order = mesh->get_q_order();

  this->addEvaluatedField(wDv);
  this->addEvaluatedField(BF);
  this->addEvaluatedField(gBF);
  this->setName("Basis Functions");
}

PHX_POST_REGISTRATION_SETUP(BasisFunctions, data, fm)
{
  this->utils.setFieldData(wDv, fm);
  this->utils.setFieldData(BF, fm);
  this->utils.setFieldData(gBF, fm);
}

PHX_EVALUATE_FIELDS(BasisFunctions, workset)
{
  apf::Mesh* m = mesh->get_apf_mesh();
  apf::FieldShape* s = mesh->get_apf_shape();

  apf::Vector3 p;
  apf::NewArray<double> bf;
  apf::NewArray<apf::Vector3> gbf;

  for (unsigned elem=0; elem < workset.size; ++elem) {
    apf::MeshEntity* e = workset.ents[elem];
    apf::MeshElement* me = apf::createMeshElement(m, e);
    for (unsigned qp=0; qp < num_qps; ++qp) {
      apf::getIntPoint(me, q_order, qp, p);
      double w = apf::getIntWeight(me, q_order, qp);
      wDv(elem, qp) = w * apf::getDV(me, p);
      apf::getBF(s, me, p, bf);
      apf::getGradBF(s, me, p, gbf);
      for (unsigned node=0; node < num_nodes; ++node) {
        BF(elem, node, qp) = bf[node];
        for (unsigned dim=0; dim < num_dims; ++dim)
          gBF(elem, node, qp, dim) = gbf[node][dim];
      }
    }
    apf::destroyMeshElement(me);
  }
}

GOAL_INSTANTIATE_ALL(BasisFunctions)

}
