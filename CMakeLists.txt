cmake_minimum_required(VERSION 2.8)
project(Goal CXX)

set(GOAL_MAJOR "0")
set(GOAL_MINOR "1")
set(GOAL_PATCH "0")

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/set_compiler_flags.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/set_options.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/find_trilinos.cmake)

add_subdirectory(src)
