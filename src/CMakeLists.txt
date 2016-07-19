configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/config.hpp")

set(HEADERS
data_types.hpp
expression.hpp
assert_param.hpp
control.hpp
dimension.hpp
layouts.hpp
traits.hpp
workset.hpp
mesh.hpp
solution_info.hpp
state_fields.hpp
mechanics.hpp
primal_problem.hpp
linear_solver.hpp
output.hpp
ev_basis_functions.hpp
ev_gather_solution.hpp
ev_dof_interpolation.hpp
ev_kinematics.hpp
ev_model_elastic.hpp
ev_model_j2.hpp
ev_first_pk.hpp
ev_mechanics_residual.hpp
ev_scatter_residual.hpp
ev_bc_dirichlet.hpp
solver.hpp
solver_continuation.hpp
solver_trapezoid.hpp
)

set(SOURCES
expression.cpp
assert_param.cpp
control.cpp
dimension.cpp
layouts.cpp
mesh.cpp
solution_info.cpp
state_fields.cpp
mechanics.cpp
mechanics_volumetric.cpp
mechanics_neumann.cpp
mechanics_dirichlet.cpp
primal_problem.cpp
linear_solver.cpp
output.cpp
ev_basis_functions.cpp
ev_gather_solution.cpp
ev_dof_interpolation.cpp
ev_kinematics.cpp
ev_model_elastic.cpp
ev_model_j2.cpp
ev_first_pk.cpp
ev_mechanics_residual.cpp
ev_scatter_residual.cpp
ev_bc_dirichlet.cpp
solver.cpp
solver_continuation.cpp
solver_trapezoid.cpp
)

add_library(goalie ${SOURCES})
target_compile_definitions(goalie PUBLIC "-DGOAL_FAD_SIZE=${GOAL_FAD_SIZE}")
target_include_directories(goalie PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
target_include_directories(goalie PUBLIC ${Trilinos_INCLUDE_DIRS})
target_include_directories(goalie PUBLIC ${Trilinos_TPL_INCLUDE_DIRS})
target_link_libraries(goalie PUBLIC scorec)
target_link_libraries(goalie PUBLIC ${Trilinos_LIBRARIES})
target_link_libraries(goalie PUBLIC ${Trilinos_TPL_LIBRARIES})
target_link_libraries(goalie PUBLIC ${Trilinos_EXTRA_LD_FLAGS})

add_executable(goal main.cpp)
target_link_libraries(goal PRIVATE goalie)