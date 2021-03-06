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

file(COPY ./meshes DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

setup_test(elast_continuation_2D)
setup_test(elast_continuation_3D)
setup_test(elast_continuation_traction_2D)
setup_test(elast_continuation_traction_3D)
setup_test(elast_continuation_mixed_2D)
setup_test(elast_continuation_mixed_3D)
setup_test(elast_continuation_temperature_3D)
setup_test(elast_continuation_body_force_2D)
setup_test(j2_continuation_2D)
setup_test(j2_continuation_3D)
setup_test(j2_continuation_mixed_2D)
setup_test(j2_continuation_mixed_3D)
setup_test(j2_continuation_mixed_traction_3D)
setup_test(j2_continuation_temperature_3D)
setup_test(j2_continuation_uniform_2D)
setup_test(j2_continuation_spr_2D)
