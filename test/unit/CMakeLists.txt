cmake_minimum_required(VERSION 3.16)

# Unit tests for MFEM-related utilities in this repo (moved to test/unit)

# Allow overriding the main shared library target name if needed
set(MAIN_SHARED_TARGET FEBioMFEMiFSI CACHE STRING "Main shared library target to link against")

include(CTest)
enable_testing()

# Helper: link a test target with FEBio libraries and set macOS rpath
# Usage: link_with_febio(<target> [FEBio::FECore FEBio::FEBioMech FEBio::FEBioXML ...])
function(link_with_febio tgt)
    if(NOT tgt)
        message(FATAL_ERROR "link_with_febio requires a <target> name")
    endif()

    # Default libraries if none provided
    if(${ARGC} GREATER 1)
        # additional libs provided; collect from ARGV1..ARGN
        set(_febio_libs "${ARGV}" )
        list(REMOVE_AT _febio_libs 0) # remove target name
    else()
        # Prefer top-level FEBio_LIBS if defined; else fall back to core libs
        if(DEFINED FEBio_LIBS)
            set(_febio_libs ${FEBio_LIBS})
        else()
            set(_febio_libs FEBio::FECore FEBio::FEBioMech FEBio::FEBioXML)
        endif()
    endif()

    target_link_libraries(${tgt} PRIVATE ${_febio_libs})

    if(APPLE)
        # Collect rpath entries for each provided FEBio lib target
        set(_rpaths)
        foreach(lib ${_febio_libs})
            list(APPEND _rpaths "$<TARGET_FILE_DIR:${lib}>")
        endforeach()
        set_target_properties(${tgt} PROPERTIES
            BUILD_WITH_INSTALL_RPATH ON
            BUILD_RPATH   "${_rpaths}"
            INSTALL_RPATH "${_rpaths}"
        )
    endif()
endfunction()

# Helper: link a test target with MFEM and all core FEBio libraries, and set macOS rpath
# Usage: link_with_mfem_and_febio(<target>)
function(link_with_mfem_and_febio tgt)
    if(NOT tgt)
        message(FATAL_ERROR "link_with_mfem_and_febio requires a <target> name")
    endif()

    # Pull libraries from top-level if available
    if(DEFINED FEBio_LIBS)
        set(_tpl_libs mfem ${FEBio_LIBS})
    else()
        set(_tpl_libs mfem FEBio::FECore FEBio::FEBioMech FEBio::FEBioXML)
    endif()

    target_link_libraries(${tgt} PRIVATE ${_tpl_libs})

    if(APPLE)
        set(_rpaths)
        foreach(lib ${_tpl_libs})
            list(APPEND _rpaths "$<TARGET_FILE_DIR:${lib}>")
        endforeach()
        set_target_properties(${tgt} PROPERTIES
            BUILD_WITH_INSTALL_RPATH ON
            BUILD_RPATH   "${_rpaths}"
            INSTALL_RPATH "${_rpaths}"
        )
    endif()
endfunction()

# Function to add a unit test executable that links MFEM and core FEBio libs
# Usage: add_tpl_unit_test(<test_name> <source_cpp_path>)
function(add_tpl_unit_test test_name src_cpp)
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
    target_include_directories("${test_name}" PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../..")

    # Link MFEM and FEBio core libraries and set rpaths on macOS
    link_with_mfem_and_febio("${test_name}")

    # Register with CTest
    add_test(NAME "${test_name}" COMMAND "${test_name}")
endfunction()

# Tests
add_tpl_unit_test(SymmetricEigenTest symmetric_eigen_test.cpp)
add_tpl_unit_test(FEBioEigenCompare febio_eigen_compare_test.cpp)
