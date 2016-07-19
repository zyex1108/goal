set(GOAL ${CMAKE_BINARY_DIR}/src/goal)

macro(copy INPUT)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/${INPUT}
    ${CMAKE_CURRENT_BINARY_DIR}/${INPUT} COPYONLY)
endmacro()

macro(setup_test TESTNAME INPUT)
  copy(${INPUT})
  add_test(${TESTNAME} ${VALGRIND} ${GOAL} ${INPUT})
endmacro()

add_subdirectory(meshes)

setup_test(elast_continuation elast_continuation.xml)
setup_test(j2_continuation j2_continuation.xml)
