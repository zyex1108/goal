set(GOALEXE ${CMAKE_BINARY_DIR}/src/goal)
set(MPIEXE mpirun)
set(MPIFLAGS -np)

macro(copy INPUT)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/${INPUT}
    ${CMAKE_CURRENT_BINARY_DIR}/${INPUT} COPYONLY)
endmacro()

macro(setup_test TESTNAME INPUT)
  copy(${INPUT})
  add_test(${TESTNAME} ${MPIEXE} ${MPIFLAGS} 4
    ${VALGRIND} ${GOALEXE} ${INPUT})
endmacro()

add_subdirectory(meshes)

setup_test(elast_continuation_2D elast_continuation_2D.xml)
setup_test(elast_continuation_3D elast_continuation_3D.xml)
setup_test(j2_continuation_2D j2_continuation_2D.xml)
setup_test(j2_continuation_3D j2_continuation_3D.xml)
setup_test(j2_continuation_2D_P2 j2_continuation_2D_P2.xml)
