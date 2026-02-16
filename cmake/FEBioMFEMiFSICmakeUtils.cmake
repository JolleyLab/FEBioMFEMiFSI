MACRO(DIRNAME result curdir)
GET_FILENAME_COMPONENT(PARENT_DIR ${curdir} DIRECTORY)

STRING(REPLACE "${PARENT_DIR}/" "" ${result} ${curdir})
ENDMACRO(DIRNAME)

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
       IF(NOT "${child}" STREQUAL "CMakeFiles")
         LIST(APPEND dirlist ${child})
       ENDIF()
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO(SUBDIRLIST)

MACRO(TARGETNAME result curdir)
GET_FILENAME_COMPONENT(PARENT_DIR ${curdir} DIRECTORY)
GET_FILENAME_COMPONENT(PARENT_PARENT_DIR ${PARENT_DIR} DIRECTORY)

STRING(REPLACE "${PARENT_PARENT_DIR}/" "" INTERMEDIATE_NAME ${curdir})
STRING(REPLACE "/" "_" ${result} ${INTERMEDIATE_NAME})
ENDMACRO(TARGETNAME)

# Function to add an executable that links to the main shared library
# Usage: ADD_iFSI_EXE(<target_name> <source_cpp_path>)
function(ADD_iFSI_EXE test_name src_cpp)
  if(NOT test_name OR NOT src_cpp)
    message(FATAL_ERROR "add_tpl_unit_test requires 2 args: <test_name> <source_cpp_path>")
  endif()

  if(NOT IS_ABSOLUTE "${src_cpp}")
    set(src_cpp "${CMAKE_CURRENT_SOURCE_DIR}/${src_cpp}")
  endif()

  if(NOT EXISTS "${src_cpp}")
    message(FATAL_ERROR "add_tpl_unit_test: source not found: ${src_cpp}")
  endif()

  add_executable("${test_name}" "${src_cpp}")

  target_compile_features("${test_name}" PUBLIC cxx_std_17)
  # Include project root for headers (test/unit is two levels deep)
  target_include_directories("${test_name}" PRIVATE "${PROJECT_SOURCE_DIR}")

  # Require the main shared target and link to it
  if(NOT DEFINED MAIN_SHARED_TARGET OR NOT TARGET ${MAIN_SHARED_TARGET})
    message(FATAL_ERROR "ADD_iFSI_EXE: MAIN_SHARED_TARGET is not defined or not built. Ensure 'FEBioMFEMiFSI' is added before '${test_name}'.")
  endif()
  target_link_libraries("${test_name}" PRIVATE ${MAIN_SHARED_TARGET})
  add_dependencies("${test_name}" ${MAIN_SHARED_TARGET})
  if(APPLE)
    # Ensure executable can find the main shared library at runtime on macOS
    set(_main_lib_rpath "$<TARGET_FILE_DIR:${MAIN_SHARED_TARGET}>")
    set_property(TARGET "${test_name}" APPEND PROPERTY BUILD_RPATH   "${_main_lib_rpath}")
    set_property(TARGET "${test_name}" APPEND PROPERTY INSTALL_RPATH "${_main_lib_rpath}")
  endif()
endfunction()