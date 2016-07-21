set(CMAKE_CXX_FLAGS "-std=c++11 -g -O2")

function(generate MODEL DIM)
  add_executable(${MODEL} ${MODEL}.cpp)
  target_include_directories(${MODEL} PRIVATE ${Gmodel_INCLUDE_DIRS})
  target_link_libraries(${MODEL} gmodel)
  add_custom_command(
    OUTPUT ${MODEL}.geo ${MODEL}.dmg
    COMMAND ${MODEL}
    DEPENDS ${MODEL}
    COMMENT "generating ${MODEL}.geo and ${MODEL}.dmg")
  if("${DIM}" EQUAL "2")
    set(GMSH_CMD ${GMSH} -2 -o ${MODEL}.msh ${MODEL}.geo)
  else()
    set(GMSH_CMD ${GMSH} -3 -optimize_netgen -o ${MODEL}.msh ${MODEL}.geo)
  endif()
  add_custom_command(
    OUTPUT ${MODEL}.msh
    COMMAND ${GMSH_CMD}
    DEPENDS ${MODEL}.geo
    COMMENT "generating ${MODEL}.msh")
  add_custom_command(
    OUTPUT ${MODEL}-serial.smb
    COMMAND ${FROMGMSH} ${MODEL}.dmg ${MODEL}.msh ${MODEL}-serial.smb
    DEPENDS ${MODEL}.dmg ${MODEL}.msh
    COMMENT "generating ${MODEL}-serial.smb")
  add_custom_command(
    OUTPUT ${MODEL}.smb
    COMMAND mpirun -np 4 ${SPLIT} ${MODEL}.dmg ${MODEL}-serial.smb ${MODEL}.smb 4
    DEPENDS ${MODEL}.dmg ${MODEL}-serial.smb
    COMMENT "generating ${MODEL}.smb")
endfunction()

macro(add_mesh MODEL DIM)
  generate(${MODEL} ${DIM})
  set(EX_MESHES ${EX_MESHES} ${MODEL}.smb)
  copy(${MODEL}.txt)
endmacro()

add_mesh(square 2)
add_mesh(plate 2)
add_mesh(cube 3)

add_custom_target(meshes
  COMMAND echo "building meshes"
  DEPENDS ${EX_MESHES})
