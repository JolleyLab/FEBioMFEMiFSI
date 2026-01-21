# CMake Housekeeping
cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_FILES_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# check to ensure apple architecture is x86_64
if(APPLE)
   execute_process(COMMAND uname -m OUTPUT_VARIABLE ARCHITECTURE OUTPUT_STRIP_TRAILING_WHITESPACE)
   if(ARCHITECTURE STREQUAL "x86_64")
      message(STATUS "Apple architecture is x86_64")
   else()
      message(FATAL_ERROR "Unsupported Apple architecture for this Plugin: ${ARCHITECTURE}")
   endif()
   set(CMAKE_C_COMPILER clang)
   set(CMAKE_CXX_COMPILER clang++)
endif()

# Set a default build type if none was specified
set(default_build_type "Release")
 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Name of your plugin
project(FEBioMFEMiFSI)

##### Find FEBio #####
message(STATUS "Looking for FEBio ...")
set(FEBio_DIR "" CACHE PATH "Path to FEBio installation.")

# Auto-detect FEBio_SOURCE_DIR based on FEBio_DIR if not explicitly set
if(NOT FEBio_SOURCE_DIR AND FEBio_DIR)
   set(FEBio_SOURCE_DIR "${FEBio_DIR}" CACHE PATH "Path to the FEBio source prefix.")
endif()

if(FEBio_DIR)
   find_package(FEBio CONFIG REQUIRED PATHS "${FEBio_DIR}" NO_DEFAULT_PATH)
else()
   find_package(FEBio CONFIG REQUIRED)
endif()

if (FEBio_FOUND)
   message(STATUS "Found FEBio in: ${FEBio_DIR}")
else()
   message(FATAL_ERROR "Could not find FEBio. Please set FEBio_DIR to the FEBio installation path.")
endif()
set(FEBio_LIBS FEBio::FECore FEBio::FEBioMech FEBio::FEBioXML)

##### Find MFEM #####
# Import MFEM. The following variables can be used to help CMake find MFEM:
#  * MFEM_DIR - absolute path to the MFEM build or install prefix.
message(STATUS "Looking for mfem ...")
if(NOT DEFINED MFEM_DIR OR MFEM_DIR STREQUAL "")
   message(STATUS "MFEM_DIR not defined, please set ...")
endif()
set(MFEM_SOURCE_DIR "" CACHE PATH "Path to the MFEM source prefix.")
# Set MFEM_SOURCE_DIR to the MFEM directory (Software/mfem) for includes
if(NOT MFEM_SOURCE_DIR)
   set(MFEM_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Software/mfem" CACHE PATH "Path to the MFEM source prefix.")
endif()
# Note: make sure to install mfem, then below will automatically find the library + header files. Otherwise set MFEM_DIR above
if (MFEM_DIR)
   find_package(MFEM CONFIG REQUIRED PATHS "${MFEM_DIR}" NO_DEFAULT_PATH)
else()
   find_package(MFEM REQUIRED NAMES MFEM)
endif()

if (MFEM_FOUND)
   message(STATUS "Found MFEM config in: ${MFEM_DIR} (version ${MFEM_VERSION})")
else()
   message(FATAL_ERROR "Could not find MFEM. Set MFEM_DIR to the folder containing MFEMConfig.cmake or add its install prefix to CMAKE_PREFIX_PATH.")
endif()
list(APPEND MFEM_SOURCE_DIR ${MFEM_INCLUDE_DIRS})

##### Set definitions and compile options #####
if(WIN32)
   add_definitions(-DWIN32 -DFECORE_DLL /MP /openmp)
elseif(APPLE)
   add_definitions(-D__APPLE__)
   #  set(CMAKE_OSX_DEPLOYMENT_TARGET 10.13)
   set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
else()
	add_definitions(-DLINUX)
   add_compile_options(-static-libstdc++ -static-libgcc -w -Wall)
    
   set(CMAKE_BUILD_RPATH_USE_LINK_PATH FALSE)
   set(CMAKE_BUILD_RPATH $ORIGIN/../lib/)
endif()

# Options file
CONFIGURE_FILE(FEBioMFEMiFSIOptions.cmake.in FEBioMFEMiFSIOptions.cmake @ONLY)

# Set include directory
include_directories(${FEBio_SOURCE_DIR} ${MFEM_SOURCE_DIR})

# Add library (include sources from subfolders)
file(GLOB_RECURSE HDR_FEBioMFEMiFSI CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
)
file(GLOB_RECURSE SRC_FEBioMFEMiFSI CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

# Exclude TPLs
list(FILTER SRC_FEBioMFEMiFSI EXCLUDE REGEX "^${CMAKE_CURRENT_SOURCE_DIR}/Software/.*")
list(FILTER HDR_FEBioMFEMiFSI EXCLUDE REGEX "^${CMAKE_CURRENT_SOURCE_DIR}/Software/.*")

# Exclude sources from the MFEMTools subdirectory; built separately
list(FILTER SRC_FEBioMFEMiFSI EXCLUDE REGEX "^${CMAKE_CURRENT_SOURCE_DIR}/MFEMTools/.*")

# Organize files by folder in IDEs
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES
    ${HDR_FEBioMFEMiFSI} ${SRC_FEBioMFEMiFSI}
)

add_library(FEBioMFEMiFSI SHARED ${HDR_FEBioMFEMiFSI} ${SRC_FEBioMFEMiFSI})
set_property(TARGET FEBioMFEMiFSI PROPERTY AUTOGEN_BUILD_DIR ${CMAKE_BINARY_DIR}/CMakeFiles/AutoGen/FEBioMFEMiFSI_autogen)

# Link FEBio Libraries (now actual paths)
target_link_libraries(FEBioMFEMiFSI PRIVATE ${FEBio_LIBS})

# Link MFEM
target_link_libraries(FEBioMFEMiFSI PRIVATE mfem)

# Propagate MFEM's compile definitions (MFEM_USE_MPI, MFEM_USE_GSLIB, etc.)
target_compile_definitions(FEBioMFEMiFSI PRIVATE $<TARGET_PROPERTY:mfem,INTERFACE_COMPILE_DEFINITIONS>)

# macOS: ensure loader can find FEBio libs at runtime
if (APPLE)
   set_target_properties(FEBioMFEMiFSI PROPERTIES
   BUILD_WITH_INSTALL_RPATH ON
   BUILD_RPATH   "${FEBio_LIB_DIR}"
   INSTALL_RPATH "${FEBio_LIB_DIR}"
   )
endif()

# Output the build type
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Option to build the MFEMTools subdirectory (off by default)
option(BUILD_MFEMTOOLS "Build MFEMTools helper/library targets" OFF)

# Optionally build tools subdirectory
if(BUILD_MFEMTOOLS)
    add_subdirectory(MFEMTools)
endif()