set(GOAL_FAD_SIZE "30" CACHE STRING "sacado derivative array size")
option(GOAL_DISABLE_CHECKS "disable sanity checks for speed" OFF)
option(GOAL_DISABLE_OPTIMIZATION "disable compiler optimization" OFF)
add_definitions(-DGOAL_FAD_SIZE=${GOAL_FAD_SIZE})
