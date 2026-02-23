ADD_iFSI_EXE(TEST_PulsatilePipeFlow PulsatilePipeFlow.cpp)

# Copy all .msh files to the build directory
file(GLOB MSH_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.msh")
foreach(MSH_FILE ${MSH_FILES})
  add_custom_command(
    TARGET TEST_PulsatilePipeFlow POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${MSH_FILE}
    ${CMAKE_CURRENT_BINARY_DIR}
  )
endforeach()