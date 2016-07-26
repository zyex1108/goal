#include "ev_elem_size.hpp"
#include "traits.hpp"
#include "layouts.hpp"
#include "workset.hpp"
#include "phx_utils.hpp"

namespace goal {

PHX_EVALUATOR_CTOR(ElemSize, p) :
  dl            (p.get<RCP<Layouts> >("Layouts")),
  size          (p.get<std::string>("Size Name"), dl->qp_scalar)
{
  num_nodes = dl->node_qp_scalar->dimension(1);
  num_qps = dl->node_qp_vector->dimension(2);
  num_dims = dl->node_qp_vector->dimension(3);

  std::string bf_name = p.get<std::string>("BF Name");
  get_grad_field(bf_name, dl, gBF);

  this->addDependentField(gBF);
  this->addEvaluatedField(size);
  this->setName("Elem Size");
}

PHX_POST_REGISTRATION_SETUP(ElemSize, data, fm)
{
  this->utils.setFieldData(gBF, fm);
  this->utils.setFieldData(size, fm);
}

PHX_EVALUATE_FIELDS(ElemSize, workset)
{
  double h=0.0;
  for (unsigned elem=0; elem < workset.size; ++elem) {
    for (unsigned qp=0; qp < num_qps; ++qp) {
      h = 0.0;
      for (unsigned i=0; i < num_dims; ++i) {
        for (unsigned node=0; node < num_nodes; ++node) {
          h += std::abs(gBF(elem,node,qp,i)/std::sqrt(num_dims));
        }
      }
      size(elem,qp) = 2.0/h;
    }
  }
}

GOAL_INSTANTIATE_ALL(ElemSize)

}
