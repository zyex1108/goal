configure_file(config_in.hpp config.hpp)

include_directories(
${Trilinos_INCLUDE_DIRS}
${Trilinos_TPL_INCLUDE_DIRS}
${CMAKE_CURRENT_BINARY_DIR}
)

set(HEADERS
expression.hpp
control.hpp
)

set(SOURCES
expression.cpp
control.cpp
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
