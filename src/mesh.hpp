#ifndef goal_mesh_hpp
#define goal_mesh_hpp

#include "data_types.hpp"

#include <apfDynamicArray.h>
#include <Teuchos_RCP.hpp>

namespace Teuchos {
class ParameterList;
}

namespace apf {
class Mesh2;
class MeshEntity;
class FieldShape;
template <class T> class NumberingOf;
typedef NumberingOf<long> GlobalNumbering;
struct StkModels;
struct Node;
}

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class Mesh
{
  public:

    Mesh(RCP<const ParameterList> p);
    ~Mesh();

    void set_num_eqs(unsigned neq) {num_eqs = neq;}
    unsigned get_num_eqs() {return num_eqs;}
    unsigned get_num_dims() {return num_dims;}
    unsigned get_ws_size() {return ws_size;}
    unsigned get_p_order() {return p_order;}
    unsigned get_q_order() {return q_order;}

    RCP<const Map> get_owned_map() {return owned_map;}
    RCP<const Map> get_overlap_map() {return overlap_map;}
    RCP<const Graph> get_owned_graph() {return owned_graph;}
    RCP<const Graph> get_overlap_graph() {return overlap_graph;}

    apf::Mesh2* get_apf_mesh() {return mesh;}
    apf::FieldShape* get_apf_shape() {return shape;}
    apf::GlobalNumbering* get_apf_numbering() {return numbering;}
    apf::DynamicArray<apf::Node> const& get_apf_nodes() {return nodes;}

    unsigned get_num_elem_qps() const;
    unsigned get_num_elem_nodes() const;
    unsigned get_num_elem_dofs() const;

    unsigned get_num_elem_sets() const;
    unsigned get_num_facet_sets() const;
    unsigned get_num_node_sets() const;

    std::string const& get_elem_set_name(const unsigned i) const;
    std::string const& get_facet_set_name(const unsigned i) const;
    std::string const& get_node_set_name(const unsigned i) const;

    unsigned get_num_worksets(const unsigned elem_set_idx);

    std::vector<apf::MeshEntity*> const& get_elems(
        std::string const& elem_set, const unsigned ws_idx);
    std::vector<apf::MeshEntity*> const& get_facets(
        std::string const& facet_set);
    std::vector<apf::Node*> const& get_nodes(
        std::string const& node_set);

    LO get_lid(apf::MeshEntity* e, const unsigned n, const unsigned eq);
    LO get_lid(apf::Node* n, const unsigned eq);

    double get_mesh_size(apf::MeshEntity* e);

    void change_p(int add);
    void update();

  private:

    RCP<const ParameterList> params;

    unsigned num_dims;
    unsigned ws_size;
    unsigned p_order;
    unsigned q_order;

    unsigned num_eqs;

    apf::Mesh2* mesh;
    apf::StkModels* sets;
    apf::FieldShape* shape;
    apf::GlobalNumbering* numbering;
    apf::DynamicArray<apf::Node> nodes;
    int elem_type;

    RCP<const Comm> comm;
    RCP<const Map> owned_map;
    RCP<const Map> overlap_map;
    RCP<Graph> owned_graph;
    RCP<Graph> overlap_graph;

    std::map<std::string, std::vector<std::vector<apf::MeshEntity*> > > elem_sets;
    std::map<std::string, std::vector<apf::MeshEntity*> > facet_sets;
    std::map<std::string, std::vector<apf::Node*> > node_sets;

    void compute_owned_map();
    void compute_overlap_map();
    void compute_graphs();

    void compute_elem_sets();
    void compute_facet_sets();
    void compute_node_sets();

};

RCP<Mesh> mesh_create(RCP<const ParameterList> p);

}

#endif
