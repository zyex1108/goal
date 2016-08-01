#include "solution_info.hpp"
#include "mesh.hpp"
#include "control.hpp"

namespace goal {

void SolutionInfo::resize(
    RCP<Mesh> mesh,
    bool enable_dynamics)
{
  double t0 = time();
  int num_vectors = 1;
  if (enable_dynamics) num_vectors = 3;
  RCP<const Map> m = mesh->get_owned_map();
  RCP<const Map> om = mesh->get_overlap_map();
  RCP<const Graph> g = mesh->get_owned_graph();
  RCP<const Graph> og = mesh->get_overlap_graph();
  exporter = rcp(new Export(om,m));
  importer = rcp(new Import(m,om));
  owned_solution = rcp(new MultiVector(m, num_vectors));
  owned_residual = rcp(new Vector(m));
  owned_jacobian = rcp(new Matrix(g));
  ovlp_solution = rcp(new MultiVector(om, num_vectors));
  ovlp_residual = rcp(new Vector(om));
  ovlp_jacobian = rcp(new Matrix(og));
  double t1 = time();
  print("solution containers resized in %f seconds", t1-t0);
}

void SolutionInfo::create_dual_vectors(
    RCP<Mesh> mesh)
{
  RCP<const Map> m = mesh->get_owned_map();
  RCP<const Map> om = mesh->get_overlap_map();
  owned_qoi = rcp(new Vector(m));
  owned_dual = rcp(new Vector(m));
  ovlp_qoi = rcp(new Vector(om));
  ovlp_dual = rcp(new Vector(om));
}

void SolutionInfo::destroy_dual_vectors()
{
  owned_qoi = Teuchos::null;
  owned_dual = Teuchos::null;
  ovlp_qoi = Teuchos::null;
  ovlp_dual = Teuchos::null;
}

void SolutionInfo::gather_solution()
{
  owned_solution->doExport(*ovlp_solution, *exporter, Tpetra::INSERT);
}

void SolutionInfo::scatter_solution()
{
  ovlp_solution->doImport(*owned_solution, *importer, Tpetra::INSERT);
}

void SolutionInfo::scatter_dual()
{
  ovlp_dual->doImport(*owned_dual, *importer, Tpetra::INSERT);
}

void SolutionInfo::gather_residual()
{
  owned_residual->doExport(*ovlp_residual, *exporter, Tpetra::ADD);
}

void SolutionInfo::gather_qoi()
{
  owned_qoi->doExport(*ovlp_qoi, *exporter, Tpetra::ADD);
}

void SolutionInfo::gather_jacobian()
{
  owned_jacobian->doExport(*ovlp_jacobian, *exporter, Tpetra::ADD);
}

RCP<SolutionInfo> sol_info_create(RCP<Mesh> m, bool enable_dynamics)
{
  RCP<SolutionInfo> s = rcp(new SolutionInfo);
  s->resize(m, enable_dynamics);
  return s;
}

}
