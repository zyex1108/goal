#include "mesh.hpp"
#include "control.hpp"
#include "assert_param.hpp"

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <apfAlbany.h>
#include <apfNumbering.h>
#include <gmi_mesh.h>

namespace goal {

static RCP<ParameterList> get_valid_params()
{
  RCP<ParameterList> p = rcp(new ParameterList);
  p->set<std::string>("geom file", "");
  p->set<std::string>("mesh file", "");
  p->set<std::string>("assoc file", "");
  p->set<unsigned>("p order", 1);
  p->set<unsigned>("q order", 1);
  p->set<unsigned>("ws size", 0);
  return p;
}

static void validate_params(RCP<const ParameterList> p)
{
  assert_param(p, "geom file");
  assert_param(p, "mesh file");
  assert_param(p, "assoc file");
  assert_param(p, "p order");
  assert_param(p, "q order");
  assert_param(p, "ws size");
  p->validateParameters(*get_valid_params(), 0);
}

static void load_mesh_from_file(
    apf::Mesh2** mesh,
    RCP<const ParameterList> p)
{
  gmi_register_mesh();
  std::string const& geom_file = p->get<std::string>("geom file");
  std::string const& mesh_file = p->get<std::string>("mesh file");
  const char* g = geom_file.c_str();
  const char* m = mesh_file.c_str();
  *mesh = apf::loadMdsMesh(g, m);
  apf::reorderMdsMesh(*mesh);
  (*mesh)->verify();
}

static apf::StkModels* read_sets(apf::Mesh* m, RCP<const ParameterList> p)
{
  apf::StkModels* sets = new apf::StkModels;
  std::string const& fn = p->get<std::string>("assoc file");
  const char* filename = fn.c_str();
  print("reading association file: %s", filename);
  static std::string const setNames[3] = {
    "node set",
    "facet set",
    "element set"};
  int d = m->getDimension();
  int dims[3] = {0, d-1, d};
  std::ifstream f(filename);
  if (!f.good())
    fail("cannot open file: %s", filename);
  std::string sline;
  int lc = 0;
  while(std::getline(f, sline)) {
    if (!sline.length())
      break;
    ++lc;
    int sdi = -1;
    for (int di=0; di < 3; ++di)
      if (sline.compare(0, setNames[di].length(), setNames[di]) == 0)
        sdi = di;
    if (sdi == -1)
      fail("invalid association line # %d:\n\t%s", lc, sline.c_str());
    int sd = dims[sdi];
    std::stringstream strs(sline.substr(setNames[sdi].length()));
    apf::StkModel* set = new apf::StkModel();
    strs >> set->stkName;
    int nents;
    strs >> nents;
    if (!strs)
      fail("invalid association line # %d:\n\t%s", lc, sline.c_str());
    for (int ei=0; ei < nents; ++ei) {
      std::string eline;
      std::getline(f, eline);
      if (!f || !eline.length())
        fail("invalid association after line # %d", lc);
      ++lc;
      std::stringstream strs2(eline);
      int mdim, mtag;
      strs2 >> mdim >> mtag;
      if (!strs2)
        fail("bad associations line # %d:\n\t%s", lc, eline.c_str());
      set->ents.push_back(m->findModelEntity(mdim, mtag));
      if (!set->ents.back())
        fail("no model entity with dim: %d and tag: %d", mdim, mtag);
    }
    sets->models[sd].push_back(set);
  }
  sets->computeInverse();
  return sets;
}

static unsigned get_elem_type(apf::Mesh* m)
{
  apf::MeshEntity* elem;
  apf::MeshIterator* it = m->begin(m->getDimension());
  elem = m->iterate(it);
  apf::Mesh::Type type = m->getType(elem);
  m->end(it);
  return type;
}

Mesh::Mesh(RCP<const ParameterList> p) :
  params(p),
  num_eqs(0),
  mesh(0),
  shape(0),
  numbering(0)
{
  validate_params(params);
  load_mesh_from_file(&mesh, params);
  sets = read_sets(mesh, params);
  elem_type = get_elem_type(mesh);
  num_dims = mesh->getDimension();
  ws_size = params->get<unsigned>("ws size");
  p_order = params->get<unsigned>("p order");
  q_order = params->get<unsigned>("q order");
  shape = apf::getHierarchic(p_order);
  comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  print(" num element sets %u", get_num_elem_sets());
  print(" num facet sets %u", get_num_facet_sets());
  print(" num node sets %u", get_num_node_sets());
}

Mesh::~Mesh()
{
  if (numbering) apf::destroyGlobalNumbering(numbering);
  if (mesh) {mesh->destroyNative(); apf::destroyMesh(mesh);}
}

unsigned Mesh::get_num_elem_qps() const
{
  return apf::countGaussPoints(elem_type,q_order);
}

unsigned Mesh::get_num_elem_nodes() const
{
  return shape->getEntityShape(elem_type)->countNodes();
}

unsigned Mesh::get_num_elem_dofs() const
{
  return get_num_elem_nodes() * num_eqs;
}

unsigned Mesh::get_num_elem_sets() const
{
  return sets->models[num_dims].size();
}

unsigned Mesh::get_num_facet_sets() const
{
  return sets->models[num_dims-1].size();
}

unsigned Mesh::get_num_node_sets() const
{
  return sets->models[0].size();
}

std::string const& Mesh::get_elem_set_name(const unsigned i) const
{
  CHECK(i < get_num_elem_sets());
  return sets->models[num_dims][i]->stkName;
}

std::string const& Mesh::get_facet_set_name(const unsigned i) const
{
  CHECK(i < get_num_facet_sets());
  return sets->models[num_dims-1][i]->stkName;
}

std::string const& Mesh::get_node_set_name(const unsigned i) const
{
  CHECK(i < get_num_node_sets());
  return sets->models[0][i]->stkName;
}

unsigned Mesh::get_num_worksets(const unsigned set_idx)
{
  std::string const& set = get_elem_set_name(set_idx);
  CHECK(elem_sets.count(set));
  return elem_sets[set].size();
}

std::vector<apf::MeshEntity*> const& Mesh::get_elems(
    std::string const& elem_set_idx, const unsigned ws_idx)
{
  return elem_sets[elem_set_idx][ws_idx];
}

std::vector<apf::MeshEntity*> const& Mesh::get_facets(
    std::string const& facet_set_idx)
{
  return facet_sets[facet_set_idx];
}

std::vector<apf::Node*> const& Mesh::get_nodes(
    std::string const& node_set_idx)
{
  return node_sets[node_set_idx];
}

static GO get_dof(const GO node, const unsigned eq, const unsigned neq)
{
  return node*neq + eq;
}

/* prevent unneeded reallocation? */
apf::NewArray<long> gids;

LO Mesh::get_lid(apf::MeshEntity* e, const unsigned n, const unsigned eq)
{
  apf::getElementNumbers(numbering, e, gids);
  GO dof = get_dof(gids[n], eq, num_eqs);
  LO lid = overlap_map->getLocalElement(dof);
  CHECK(lid >=0);
  return lid;
}

LO Mesh::get_lid(apf::Node* node, const unsigned eq)
{
  CHECK(mesh->isOwned(node->entity));
  long n = apf::getNumber(numbering, *node);
  GO dof = get_dof(n, eq, num_eqs);
  LO lid = owned_map->getLocalElement(dof);
  CHECK(lid >= 0);
  return lid;
}

void Mesh::compute_owned_map()
{
  if (numbering) apf::destroyGlobalNumbering(numbering);
  numbering = apf::makeGlobal(apf::numberOwnedNodes(mesh,"n",shape));
  apf::DynamicArray<apf::Node> owned;
  apf::getNodes(numbering, owned);
  unsigned num_owned_nodes = owned.getSize();
  Teuchos::Array<GO> indices(num_eqs*num_owned_nodes);
  for (unsigned i=0; i < num_owned_nodes; ++i) {
    GO gid = apf::getNumber(numbering, owned[i]);
    for (unsigned j=0; j < num_eqs; ++j)
      indices[get_dof(i,j,num_eqs)] = get_dof(gid,j,num_eqs);
  }
  owned_map = Tpetra::createNonContigMap<LO,GO>(indices,comm);
  apf::synchronize(numbering);
}

void Mesh::compute_overlap_map()
{
  apf::Numbering* overlap = apf::numberOverlapNodes(mesh,"o",shape);
  apf::getNodes(overlap,nodes);
  unsigned num_overlap_nodes = nodes.getSize();
  Teuchos::Array<GO> indices(num_eqs*num_overlap_nodes);
  for (unsigned i=0; i < num_overlap_nodes; ++i) {
    GO gid = apf::getNumber(numbering,nodes[i]);
    for (unsigned j=0; j < num_eqs; ++j)
      indices[get_dof(i,j,num_eqs)] = get_dof(gid,j,num_eqs);
  }
  overlap_map = Tpetra::createNonContigMap<LO,GO>(indices,comm);
  apf::destroyNumbering(overlap);
}

static unsigned estimate_bdwth(const unsigned neqs, const unsigned ndims)
{
  unsigned est_bdwth = 0;
  switch (ndims) {
    case 0: est_bdwth = 1*neqs; break;
    case 1: est_bdwth = 3*neqs; break;
    case 2: est_bdwth = 9*neqs; break;
    case 3: est_bdwth = 27*neqs; break;
    default: fail("bad number of dimensions: %d", ndims);
  }
  return est_bdwth;
}

void Mesh::compute_graphs()
{
  overlap_graph = rcp(new Graph(overlap_map, get_num_elem_dofs()));
  apf::MeshEntity* elem;
  apf::MeshIterator* elems = mesh->begin(num_dims);
  while ((elem = mesh->iterate(elems))) {
    apf::NewArray<long> cell_nodes;
    unsigned nnodes = apf::getElementNumbers(numbering,elem,cell_nodes);
    for (unsigned i=0; i < nnodes; ++i) {
      for (unsigned j=0; j < num_eqs; ++j) {
        GO row = get_dof(cell_nodes[i],j,num_eqs);
        for (unsigned k=0; k < nnodes; ++k) {
          for (unsigned l=0; l < num_eqs; ++l) {
            GO col = get_dof(cell_nodes[k],l,num_eqs);
            Teuchos::ArrayView<GO> colav = Teuchos::arrayView(&col,1);
            overlap_graph->insertGlobalIndices(row,colav);
  }}}}}
  mesh->end(elems);
  overlap_graph->fillComplete();
  unsigned r = estimate_bdwth(num_eqs,num_dims);
  owned_graph = rcp(new Graph(owned_map,r));
  RCP<Export> exporter = rcp(new Export(overlap_map,owned_map));
  owned_graph->doExport(*overlap_graph,*exporter,Tpetra::INSERT);
  owned_graph->fillComplete();
}

void Mesh::compute_elem_sets()
{
  unsigned nes = get_num_elem_sets();
  for (unsigned i=0; i < nes; ++i)
    elem_sets[get_elem_set_name(i)].resize(0);
  apf::MeshEntity* elem;
  apf::MeshIterator* it = mesh->begin(num_dims);
  std::map<std::string, std::vector<apf::MeshEntity*> > map;
  for (unsigned i=0; i < nes; ++i)
    map[get_elem_set_name(i)].resize(0);
  while ((elem = mesh->iterate(it))) {
    apf::ModelEntity* mr = mesh->toModel(elem);
    apf::StkModel* stkm = sets->invMaps[num_dims][mr];
    std::string const& name = stkm->stkName;
    if (map[name].size() >= ws_size) {
      elem_sets[name].push_back(map[name]);
      map[name].clear();
    }
    map[name].push_back(elem);
  }
  mesh->end(it);
  for (unsigned i=0; i < nes; ++i) {
    std::string const& name = get_elem_set_name(i);
    elem_sets[name].push_back(map[name]);
  }
}

void Mesh::compute_facet_sets()
{
  unsigned nfs = sets->models[num_dims-1].size();
  for (unsigned i=0; i < nfs; ++i)
    facet_sets[get_facet_set_name(i)].resize(0);
  apf::MeshIterator* it = mesh->begin(num_dims-1);
  apf::MeshEntity* facet;
  while ((facet = mesh->iterate(it))) {
    apf::ModelEntity* me = mesh->toModel(facet);
    if (!sets->invMaps[num_dims-1].count(me))
      continue;
    apf::StkModel* fs = sets->invMaps[num_dims-1][me];
    std::string const& fsn = fs->stkName;
    apf::Up adjElems;
    mesh->getUp(facet, adjElems);
    CHECK(adjElems.n == 1);
    facet_sets[fsn].push_back(facet);
  }
  mesh->end(it);
}

void Mesh::compute_node_sets()
{
  unsigned nds = sets->models[0].size();
  for (unsigned i=0; i < nds; ++i)
    node_sets[get_node_set_name(i)].resize(0);
  for (unsigned i=0; i < nodes.size(); ++i) {
    apf::Node* node = &(nodes[i]);
    apf::MeshEntity* e = node->entity;
    if (!mesh->isOwned(e))
      continue;
    std::set<apf::StkModel*> mset;
    apf::collectEntityModels(
        mesh, sets->invMaps[0], mesh->toModel(e), mset);
    if (mset.empty())
      continue;
    APF_ITERATE(std::set<apf::StkModel*>, mset, mit) {
      apf::StkModel* ds = *mit;
      std::string const& dsn = ds->stkName;
      node_sets[dsn].push_back(node);
    }
  }
}

void Mesh::update()
{
  double t0 = time();
  compute_owned_map();
  compute_overlap_map();
  compute_graphs();
  compute_elem_sets();
  compute_facet_sets();
  compute_node_sets();
  double t1 = time();
  print("mesh updated in %f seconds", t1-t0);
}

RCP<Mesh> mesh_create(RCP<const ParameterList> p)
{
  RCP<const ParameterList> mp = rcpFromRef(p->sublist("mesh"));
  return rcp(new Mesh(mp));
}

}
