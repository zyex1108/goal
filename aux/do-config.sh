# Modify these paths for your system
SCORECDIR=
TRIDIR=

cmake \
 -D CMAKE_CXX_COMPILER="mpicxx" \
 -D SCOREC_PREFIX=$SCORECDIR \
 -D Trilinos_PREFIX=$TRIDIR \
 -D GOAL_FAD_SIZE=16 \
 -D GOAL_TESTING=ON \
 -D GOAL_OPTIMIZE=ON \
 -D GOAL_SYMBOLS=ON \
 -D GOAL_DISABLE_CHECKS=OFF \
..
