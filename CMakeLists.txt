cmake_minimum_required(VERSION 2.8)
project(Goal LANGUAGES CXX)

set(GOAL_MAJOR "0")
set(GOAL_MINOR "17")
set(GOAL_PATCH "0")

set(SCOREC_PREFIX "" CACHE PATH "Scorec installation directory")
set(Trilinos_PREFIX "" CACHE PATH "Trilinos installation directory")
set(Valgrind_PREFIX "" CACHE PATH "Valgrind installation directory")
set(GOAL_FAD_SIZE "16" CACHE STRING "Maximum Sacado derivative array size")
set(GOAL_TESTING OFF CACHE STRING "Build and enable tests")
set(GOAL_VALGRIND OFF CACHE STRING "Run tests under valgrind")
set(GOAL_EXAMPLES OFF CACHE STRING "Clone and enable the examples repository")
option(GOAL_OPTIMIZE "Compile with optimizations" ON)
option(GOAL_SYMBOLS "Compile with symbols" ON)
option(GOAL_DISABLE_CHECKS "Disable basic sanity checks for speed" OFF)

message(STATUS "GOAL_FAD_SIZE: ${GOAL_FAD_SIZE}")
message(STATUS " maximum Sacado derivative array size")
message(STATUS "GOAL_OPTIMIZE: ${GOAL_OPTIMIZE}")
message(STATUS " compile with optimization flags")
message(STATUS "GOAL_SYMBOLS: ${GOAL_SYMBOLS}")
message(STATUS " compile with debug symbols")
message(STATUS "GOAL_DISABLE_CHECKS: ${GOAL_DISABLE_CHECKS}")
message(STATUS " disbale basic sanity checks for speed")
message(STATUS "GOAL_TESTING: ${GOAL_TESTING}")
message(STATUS " build and enable tests")
message(STATUS "GOAL_VALGRIND: ${GOAL_VALGRIND}")
message(STATUS " run tests under valgrind")
message(STATUS "GOAL_EXAMPLES: ${GOAL_EXAMPLES}")
message(STATUS " clone and enable the examples repository")

set(FLAGS "-std=c++11")
if(GOAL_OPTIMIZE)
  set(FLAGS "${FLAGS} -O2")
else()
  set(FLAGS "${FLAGS} -O0")
endif()
if(GOAL_SYMBOLS)
  set(FLAGS "${FLAGS} -g")
endif()
set(FLAGS "${FLAGS} -Wall -Wpedantic -Werror")
set(FLAGS "${FLAGS} -fno-omit-frame-pointer")
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(FLAGS "${FLAGS} -fmax-errors=1")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  set(FLAGS "${FLAGS} -ferror-limit=1")
endif()
set(CMAKE_CXX_FLAGS "${FLAGS}")

set(CMAKE_SKIP_BUILD_RPATH False)
set(CMAKE_BUILD_WITH_INSTALL_RPATH False)
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
  "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True)

find_package(SCOREC 2.0.0 REQUIRED PATHS ${SCOREC_PREFIX})
find_package(Trilinos 12.0.0 REQUIRED PATHS ${Trilinos_PREFIX})

list(REMOVE_DUPLICATES Trilinos_INCLUDE_DIRS)
list(REMOVE_DUPLICATES Trilinos_TPL_INCLUDE_DIRS)
list(REMOVE_DUPLICATES Trilinos_LIBRARIES)
list(REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)

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

list(FIND Trilinos_TPL_LIST MPI MPIListIdx)
if(NOT MPIListIdx GREATER -1)
message(FATAL_ERROR "Trilinos: mpi not enabled")
endif()

list(FIND Trilinos_TPL_LIST Boost BoostIdx)
if(NOT BoostIdx GREATER -1)
message(FATAL_ERROR "Trilinos: boost is not enabled")
endif()

message(STATUS "found Trilinos: ${Trilinos_DIR} (${Trilinos_VERSION})")
message(STATUS "found SCOREC: ${SCOREC_DIR} (${SCOREC_VERSION})")

add_subdirectory(src)

if(GOAL_TESTING)
  include(CTest)
  if(GOAL_VALGRIND)
    message(STATUS "found Valgrind: ${VALGRIND}")
  endif()
  add_subdirectory(tests)
endif()

if(GOAL_EXAMPLES)
  if(NOT GOAL_TESTING)
    message(FATAL_ERROR "GOAL_EXAMPLES=ON requires GOAL_TESTING=ON!")
  endif()
  if(NOT IS_DIRECTORY ${CMAKE_SOURCE_DIR}/examples)
    set(GIT_CMD "git")
    set(GIT_ARG1 "clone")
    set(GIT_ARG2 "git@github.com:bgranzow/goal-examples")
    set(GIT_ARG3 "examples")
    execute_process(
      COMMAND ${GIT_CMD} ${GIT_ARG1} ${GIT_ARG2} ${GIT_ARG3}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
endif()
