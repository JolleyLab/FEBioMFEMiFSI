cmake_minimum_required(VERSION 3.16)

# Name of the main shared library target from the top-level CMakeLists
set(MAIN_SHARED_TARGET FEBioMFEMiFSI CACHE STRING "Main shared library target to link against")

# Function to add a tool executable
# Usage: add_mfem_tool(<exe_name> <source_cpp_path>)
function(add_mfem_tool exe_name src_cpp)
    if(NOT exe_name OR NOT src_cpp)
        message(FATAL_ERROR "add_mfem_tool requires 2 args: <exe_name> <source_cpp_path>")
    endif()

    if(NOT IS_ABSOLUTE "${src_cpp}")
        set(src_cpp "${CMAKE_CURRENT_SOURCE_DIR}/${src_cpp}")
    endif()

    if(NOT EXISTS "${src_cpp}")
        message(FATAL_ERROR "add_mfem_tool: source not found: ${src_cpp}")
    endif()

    add_executable("${exe_name}" "${src_cpp}")

    target_compile_features("${exe_name}" PUBLIC cxx_std_17)
    target_include_directories("${exe_name}" PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

    # Link
    target_link_libraries("${exe_name}" PRIVATE mfem)
    target_link_libraries("${exe_name}" PRIVATE ${MAIN_SHARED_TARGET})
endfunction()

# Tools
# add_mfem_tool(<exe_name> <source_cpp_path>)
add_mfem_tool(MeshElemSize mesh_elem_size.cpp)
add_mfem_tool(CheckMeshSizeRatio check_mesh_size_ratio.cpp)