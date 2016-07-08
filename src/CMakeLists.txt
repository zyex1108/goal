configure_file(config_in.hpp config.hpp)

include_directories(
${Trilinos_INCLUDE_DIRS}
${Trilinos_TPL_INCLUDE_DIRS}
${CMAKE_CURRENT_BINARY_DIR}
)

set(HEADERS
data_types.hpp
expression.hpp
assert_param.hpp
control.hpp
dimension.hpp
layouts.hpp
traits.hpp
mesh.hpp
state_fields.hpp
mechanics.hpp
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
state_fields.cpp
mechanics.cpp
solver.cpp
solver_continuation.cpp
solver_trapezoid.cpp
)

add_library(goallib ${SOURCES} ${HEADERS})

add_executable(goal main.cpp)
target_link_libraries(
goal
goallib
${Trilinos_LIBRARIES}
${Trilinos_TPL_LIBRARIES}
${Trilinos_EXTRA_LD_FLAGS}
)
