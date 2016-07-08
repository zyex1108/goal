set(CMAKE_PREFIX_PATH ${TRILINOS_DIR} ${CMAKE_PREFIX_PATH})

find_package(Trilinos REQUIRED QUIET)

list(REMOVE_DUPLICATES Trilinos_INCLUDE_DIRS)
list(REMOVE_DUPLICATES Trilinos_TPL_INCLUDE_DIRS)
list(REMOVE_DUPLICATES Trilinos_LIBRARIES)
list(REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)

list(FIND Trilinos_TPL_LIST MPI MPIListIdx)
if(NOT MPIListIdx GREATER -1)
message(FATAL_ERROR "Trilinos: mpi not enabled")
endif()

list(FIND Trilinos_PACKAGE_LIST Pamgen PamgenListIdx)
if(NOT PamgenListIdx GREATER -1)
message(FATAL_ERROR "Trilinos: pamgen not enabled")
endif()

list(FIND Trilinos_PACKAGE_LIST Phalanx PhalanxListIdx)
if(NOT PhalanxListIdx GREATER -1)
message(FATAL_ERROR "Trilinos: phalanx not enabled")
endif()

list(FIND Trilinos_PACKAGE_LIST Belos BelosListIdx)
if(NOT BelosListIdx GREATER -1)
message(FATAL_ERROR "Trilinos: belos not enabled")
endif()

list(FIND Trilinos_PACKAGE_LIST Ifpack2 Ifpack2Idx)
if(NOT Ifpack2Idx GREATER -1)
message(FATAL_ERROR "Trilinos: ifpack2 not enabled")
endif()

list(FIND Trilinos_TPL_LIST Boost BoostIdx)
if(NOT BoostIdx GREATER -1)
message(FATAL_ERROR "Trilinos: boost is not enabled")
endif()

list(FIND Trilinos_PACKAGE_LIST SCOREC SCORECListIdx)
if(NOT SCORECListIdx GREATER -1)
message(FATAL_ERROR "Trilinos: scorec is not enabled")
endif()

list(FIND Trilinos_PACKAGE_LIST MueLu MueLuIdx)
list(FIND Trilinos_PACKAGE_LIST Amesos2 Amesos2Idx)
if((MueLuIdx GREATER -1) AND (Amesos2Idx GREATER -1))
set(GOAL_ENABLE_AMG ON)
endif()

message("-- Found Trilinos: ${Trilinos_DIR} (${Trilinos_VERSION})")
