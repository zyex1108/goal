set(GOALEXE ${CMAKE_BINARY_DIR}/src/goal)
set(MPIEXE mpirun)
set(MPIFLAGS -np)

macro(copy INPUT)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/${INPUT}
    ${CMAKE_CURRENT_BINARY_DIR}/${INPUT} COPYONLY)
endmacro()

macro(setup_test INPUT)
  set(INPUTXML "${INPUT}.xml")
  copy(${INPUTXML})
  add_test(${INPUT} ${MPIEXE} ${MPIFLAGS} 4
    ${VALGRIND} ${GOALEXE} ${INPUTXML})
endmacro()

macro(setup_serial_test INPUT)
  set(INPUTXML "${INPUT}.xml")
  copy(${INPUTXML})
  add_test(${INPUT} ${VALGRIND} ${GOALEXE} ${INPUTXML})
endmacro()

add_subdirectory(meshes)

setup_serial_test(elast_continuation_2D_serial)
setup_serial_test(elast_continuation_3D_serial)

setup_test(elast_continuation_2D)
setup_test(elast_continuation_3D)
setup_test(elast_continuation_mixed_2D)
setup_test(elast_continuation_mixed_plate)
setup_test(j2_continuation_2D)
setup_test(j2_continuation_3D)
setup_test(j2_continuation_2D_P2)
setup_test(j2_continuation_mixed_2D)
setup_test(j2_continuation_mixed_plate)
