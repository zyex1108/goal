set(CID ${CMAKE_CXX_COMPILER_ID})

set(FLAGS "-std=c++11")

if (GOAL_DISABLE_OPTIMIZATION)
set(FLAGS "${FLAGS} -g -O0")
else()
set(FLAGS "${FLAGS} -g -O2")
endif()

set(FLAGS "${FLAGS} -Werror -Wall -Wpedantic")
set(FLAGS "${FLAGS} -Wno-inconsistent-missing-override")
set(FLAGS "${FLAGS} -fno-omit-frame-pointer")

if (${CID} STREQUAL "GNU")
set(FLAGS "${FLAGS} -fmax-errors=1")
elseif(${CID} STREQUAL "Clang")
set(FLAGS "${FLAGS} -ferror-limit=1")
else()
message(FATAL_ERROR "Unsupported compiler ${CID}")
endif()

set(CMAKE_CXX_FLAGS ${FLAGS})

unset(CID)
