# Modify these paths for your system.
TRILINSTALLDIR=
TRILSRCDIR=
MPIDIR=
PARMETISDIR=
BOOSTDIR=
HAVE_LL=

cmake \
\
 -D CMAKE_INSTALL_PREFIX:PATH=$TRILINSTALLDIR \
 -D CMAKE_BUILD_TYPE:STRING=NONE \
 -D CMAKE_C_FLAGS:STRING="-g -O2" \
 -D CMAKE_CXX_FLAGS:STRING="-g -O2 -Wno-deprecated-declarations -Wno-sign-compare" \
 -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
\
 -D BUILD_SHARED_LIBS:BOOL=OFF \
\
 -D Trilinos_ENABLE_SECONDARY_TESTED_CODE:BOOL=OFF \
 -D Trilinos_ENABLE_EXPORT_MAKEFILES:BOOL=OFF \
 -D Trilinos_ASSERT_MISSING_PACKAGES:BOOL=OFF \
 -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
 -D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="" \
\
 -D Trilinos_ENABLE_Teuchos:BOOL=ON \
 -D Trilinos_ENABLE_Shards:BOOL=ON \
 -D Trilinos_ENABLE_Sacado:BOOL=ON \
 -D Trilinos_ENABLE_Belos:BOOL=ON \
 -D Trilinos_ENABLE_Intrepid2:BOOL=ON \
 -D Trilinos_ENABLE_Tpetra:BOOL=ON \
 -D Trilinos_ENABLE_Ifpack2:BOOL=ON \
\
 -D Trilinos_ENABLE_Kokkos:BOOL=ON \
 -D Trilinos_ENABLE_KokkosCore:BOOL=ON \
 -D Kokkos_ENABLE_Serial:BOOL=ON \
 -D Kokkos_ENABLE_OpenMP:BOOL=OFF \
 -D Kokkos_ENABLE_Pthread:BOOL=OFF \
\
 -D Trilinos_ENABLE_Phalanx:BOOL=ON \
 -D Phalanx_INDEX_SIZE_TYPE:STRING="INT" \
 -D Phalanx_SHOW_DEPRECATED_WARNINGS:BOOL=OFF \
 -D Phalanx_KOKKOS_DEVICE_TYPE:STRING="SERIAL" \
\
 -D Zoltan_ENABLE_ULLONG_IDS:BOOL=ON \
 -D Teuchos_ENABLE_LONG_LONG_INT:BOOL=$HAVE_LL \
\
 -D Trilinos_ENABLE_SCOREC:BOOL=ON \
 -D PCU_COMPRESS:BOOL=ON \
 -D SCOREC_DISABLE_STRONG_WARNINGS:BOOL=ON \
\
 -D Trilinos_ENABLE_Zoltan2:BOOL=OFF \
 -D Trilinos_ENABLE_MueLu:BOOL=OFF \
 -D Trilinos_ENABLE_Amesos2:BOOL=OFF \
 -D Trilinos_ENABLE_NOX:BOOL=OFF \
 -D Trilinos_ENABLE_Stratimikos:BOOL=OFF \
 -D Trilinos_ENABLE_Thyra:BOOL=OFF \
 -D Trilinos_ENABLE_Rythmos:BOOL=OFF \
 -D Trilinos_ENABLE_Stokhos:BOOL=OFF \
 -D Trilinos_ENABLE_Piro:BOOL=OFF \
 -D Trilinos_ENABLE_Teko:BOOL=OFF \
 -D Trilinos_ENABLE_STKIO:BOOL=OFF \
 -D Trilinos_ENABLE_STKMesh:BOOL=OFF \
 -D Trilinos_ENABLE_SEACAS:BOOL=OFF \
 -D Trilinos_ENABLE_Epetra:BOOL=OFF \
 -D Trilinos_ENABLE_EpetraExt:BOOL=OFF \
 -D Trilinos_ENABLE_Ifpack:BOOL=OFF \
 -D Trilinos_ENABLE_AztecOO:BOOL=OFF \
 -D Trilinos_ENABLE_Amesos:BOOL=OFF \
 -D Trilinos_ENABLE_Anasazi:BOOL=OFF \
 -D Trilinos_ENABLE_ML:BOOL=OFF \
 -D Trilinos_ENABLE_Intrepid:BOOL=OFF \
\
 -D TPL_ENABLE_MPI:BOOL=ON \
 -D MPI_BASE_DIR:PATH=$MPIDIR \
\
 -D TPL_ENABLE_Boost:BOOL=ON \
 -D Boost_INCLUDE_DIRS:FILEPATH="$BOOSTDIR/include" \
 -D Boost_LIBRARY_DIRS:FILEPATH="$BOOSTDIR/lib" \
\
 -D TPL_ENABLE_ParMETIS:STRING=ON \
 -D ParMETIS_INCLUDE_DIRS:PATH="$PARMETISDIR/include" \
 -D ParMETIS_LIBRARY_DIRS:PATH="$PARMETISDIR/lib" \
\
 -D TPL_ENABLE_METIS:STRING=ON \
 -D METIS_INCLUDE_DIRS:PATH="$PARMETISDIR/include" \
 -D METIS_LIBRARY_DIRS:PATH="$PARMETISDIR/lib" \
\
 -D TPL_ENABLE_Netcdf:BOOL=OFF \
 -D TPL_ENABLE_HDF5:BOOL=OFF \
 -D TPL_ENABLE_SuperLU:BOOL=OFF \
 -D TPL_ENABLE_X11:BOOL=OFF \
 -D TPL_ENABLE_Matio:BOOL=OFF \
\
 -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
 -D Tpetra_INST_FLOAT:BOOL=OFF \
 -D Tpetra_INST_INT_INT:BOOL=ON \
 -D Tpetra_INST_DOUBLE:BOOL=ON \
 -D Tpetra_INST_COMPLEX_FLOAT:BOOL=OFF \
 -D Tpetra_INST_COMPLEX_DOUBLE:BOOL=OFF \
 -D Tpetra_INST_INT_LONG:BOOL=OFF \
 -D Tpetra_INST_INT_UNSIGNED:BOOL=OFF \
 -D Tpetra_INST_INT_LONG_LONG:BOOL=$HAVE_LL \
\
$TRILSRCDIR