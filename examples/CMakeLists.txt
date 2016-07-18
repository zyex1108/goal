include(CTest)

set(CMAKE_CXX_FLAGS "-std=c++11 -g -O2")

function(generate MODEL_NAME MODEL_DIM)
  add_executable(gen_${MODEL_NAME} ${MODEL_NAME}.cpp)
  target_include_directories(gen_${MODEL_NAME}
    PRIVATE ${Gmodel_INCLUDE_DIRS})
  target_link_libraries(gen_${MODEL_NAME} gmodel)
  add_custom_command(
    OUTPUT ${MODEL_NAME}.geo ${MODEL_NAME}.dmg
    COMMAND ./gen_${MODEL_NAME}
    DEPENDS gen_${MODEL_NAME}
    COMMENT "generating ${MODEL_NAME}.geo and ${MODEL_NAME}.dmg")
  if("${MODEL_DIM}" EQUAL "2")
    set(GMSH_CMD ${GMSH} -2 -o ${MODEL_NAME}.msh ${MODEL_NAME}.geo)
  else()
    set(GMSH_CMD ${GMSH} -3 -optimize_netgen -o ${MODEL_NAME}.msh ${MODEL_NAME}.geo)
  endif()
  add_custom_command(
    OUTPUT ${MODEL_NAME}.msh
    COMMAND ${GMSH_CMD}
    DEPENDS ${MODEL_NAME}.geo
    COMMENT "generating ${MODEL_NAME}.msh")
  add_custom_command(
    OUTPUT ${MODEL_NAME}.smb
    COMMAND ${FROMGMSH} ${MODEL_NAME}.dmg ${MODEL_NAME}.msh ${MODEL_NAME}.smb
    DEPENDS ${MODEL_NAME}.dmg ${MODEL_NAME}.msh
    COMMENT "generating ${MODEL_NAME}.smb")
endfunction()

generate(square 2)
set(EX_MESHES square.smb)

add_custom_target(meshes
  COMMAND ${CMAKE_TEST_COMMAND}
  DEPENDS ${EX_MESHES})
