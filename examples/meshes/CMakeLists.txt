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
    OUTPUT ${MODEL}.smb
    COMMAND ${FROMGMSH} ${MODEL}.dmg ${MODEL}.msh ${MODEL}.smb
    DEPENDS ${MODEL}.dmg ${MODEL}.msh
    COMMENT "generating ${MODEL_NAME}.smb")
endfunction()

generate(square 2)
set(EX_MESHES ${EX_MESHES} square.smb)
generate(kirsch 2)
set(EX_MESHES ${EX_MESHES} kirsch.smb)

add_custom_target(meshes
  COMMAND echo "building meshes"
  DEPENDS ${EX_MESHES})
