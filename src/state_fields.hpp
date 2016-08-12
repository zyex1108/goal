#ifndef goal_state_fields_hpp
#define goal_state_fields_hpp

#include <Teuchos_RCP.hpp>
#include <Intrepid2_MiniTensor.h>

namespace apf {
class Mesh;
class Field;
class MeshEntity;
}

namespace goal {

using Teuchos::RCP;
using apf::MeshEntity;

class Mesh;

enum FieldType {SCALAR, VECTOR, TENSOR};

class StateFields
{
  public:

    StateFields(RCP<Mesh> m);
    
    void add(char const* n, FieldType t, bool save=false, bool I=false);

    template <typename T>
    void set_scalar(
        char const* name,
        apf::MeshEntity* e,
        unsigned qp,
        T const& v);

    template <typename T>
    void set_vector(
        char const* name,
        apf::MeshEntity* e,
        unsigned qp,
        Intrepid2::Vector<T> const& v);

    template <typename T>
    void set_tensor(
        char const* name,
        apf::MeshEntity* e,
        unsigned qp,
        Intrepid2::Tensor<T> const& v);

    template <typename T>
    void get_scalar(
        char const* name,
        apf::MeshEntity* e,
        unsigned qp,
        T& v);

    template <typename T>
    void get_vector(
        char const* name,
        apf::MeshEntity* e,
        unsigned qp,
        Intrepid2::Vector<T>& v);

    template <typename T>
    void get_tensor(
        char const* name,
        apf::MeshEntity* e,
        unsigned qp,
        Intrepid2::Tensor<T>& v);

    void project();

    void update();

  private:

    RCP<Mesh> mesh;
    apf::Mesh* apf_mesh;
    std::vector<apf::Field*> states;
    std::vector<apf::Field*> old_states;

};

}

#endif
