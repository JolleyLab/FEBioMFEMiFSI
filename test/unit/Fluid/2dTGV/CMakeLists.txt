ADD_iFSI_EXE(TEST_2dTGV 2dTGV.cpp)

# Copy tgv.msh to the build directory
add_custom_command(
  TARGET TEST_2dTGV POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  ${CMAKE_CURRENT_SOURCE_DIR}/tgv.msh
  ${CMAKE_CURRENT_BINARY_DIR}/tgv.msh
)