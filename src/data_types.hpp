#ifndef goal_data_types_hpp
#define goal_data_types_hpp

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Phalanx_KokkosDeviceTypes.hpp>
#include <Sacado_Fad_SLFad.hpp>

namespace goal {

typedef double ST;
typedef int LO;
typedef long long GO;
typedef Sacado::Fad::SLFad<ST, GOAL_FAD_SIZE> FadType;
typedef Teuchos::Comm<int> Comm;
typedef Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device> KNode;
typedef Tpetra::Map<LO, GO, KNode> Map;
typedef Tpetra::CrsGraph<LO, GO, KNode> Graph;
typedef Tpetra::Export<LO, GO, KNode> Export;
typedef Tpetra::Import<LO, GO, KNode> Import;
typedef Tpetra::Vector<ST, LO, GO, KNode> Vector;
typedef Tpetra::MultiVector<ST, LO, GO, KNode> MultiVector;
typedef Tpetra::CrsMatrix<ST, LO, GO, KNode> Matrix;
typedef Tpetra::MatrixMarket::Writer<Matrix> MM_Writer;

}

#endif
