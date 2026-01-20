#!/bin/bash

# FEBioMFEMiFSI Build Script
# This script builds the FEBioMFEMiFSI project. It assumes FEBio is built and 
# the path to the FEBio installation is set via the variable FEBio_DIR

set -e  # Exit on error

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Detect number of processors (cross-platform)
if command -v nproc >/dev/null 2>&1; then
    NUM_PROCS=$(nproc)
elif command -v sysctl >/dev/null 2>&1; then
    NUM_PROCS=$(sysctl -n hw.ncpu)
else
    NUM_PROCS=1
fi

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="${SCRIPT_DIR}"
SOFTWARE_DIR="${PROJECT_ROOT}/Software"
MFEM_DIR="${SOFTWARE_DIR}/mfem"
MFEM_BUILD_DIR="${MFEM_DIR}/build"
GSLIB_DIR="${SOFTWARE_DIR}/gslib"
GSLIB_BUILD_DIR="${GSLIB_DIR}/build"
HYPRE_DIR="${SOFTWARE_DIR}/hypre"
HYPRE_BUILD_DIR="${HYPRE_DIR}/src/hypre/"
METIS_DIR="${SOFTWARE_DIR}/metis"
METIS_BUILD_DIR="${METIS_DIR}"

# Compiler override placeholders; computed after parsing
CMAKE_COMPILER_ARGS=""
APPLE_CC=""
APPLE_CXX=""
OVERRIDE_CC=""
OVERRIDE_CXX=""
CC_CHOICE=""
CXX_CHOICE=""

# Function to print usage
usage() {
    echo -e "${YELLOW}Usage: $0 [OPTIONS]${NC}"
    echo "Required Options:"
    echo "  --febio-dir PATH            Path to FEBio installation directory (REQUIRED)"
    echo ""
    echo "Optional Options:"
    echo "  --build-type TYPE           Build type: Debug or Release (default: Release)"
    echo "  --gslib-cflags FLAGS        Additional C flags for GSLIB (default flags: -fPIC)"
    echo "  --hypre-cflags FLAGS        Additional C flags for HYPRE (default flags: -fPIC)"
    echo "  --hypre-cxxflags FLAGS      Additional C++ flags for HYPRE (default flags: -fPIC)"
    echo "  --metis-cflags FLAGS        Additional C flags for METIS (default flags: -fPIC)"
    echo "  --metis-cxxflags FLAGS      Additional C++ flags for METIS (default flags: -fPIC)"
    echo "  --mfem-cxxflags FLAGS       Additional C++ flags for MFEM (default flags: -fPIC)"
    echo "  --cc CC                     Override C compiler for CMake builds (default: platform-specific)"
    echo "  --cxx CXX                   Override C++ compiler for CMake builds (default: platform-specific)"
    echo "  --mfem-options STR          Pass additional options to MFEM CMake (e.g., '-DMFEM_USE_MPI=ON')"
    echo "  --update-ifsi               Update iFSI subtree from remote (requires clean working tree)"
    echo "  --clean                     Clean and rebuild both MFEM and FEBioMFEMiFSI plugin"
    echo "  --clean-mfem                Clean and rebuild MFEM only"
    echo "  --clean-plugin              Clean and rebuild FEBioMFEMiFSI plugin only"
    echo "  --help                      Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 --febio-dir /path/to/FEBio --febio-source-dir /path/to/FEBio/source --gslib-cflags '-O3' --mfem-cxxflags '-O3 -march=native'"
    exit 0
}

# Parse command line arguments
FEBio_DIR=""
MFEM_OPTIONS=""
GSLIB_CFLAGS=""
HYPRE_CFLAGS=""
HYPRE_CXXFLAGS=""
METIS_CFLAGS=""
METIS_CXXFLAGS=""
MFEM_CXXFLAGS=""
CLEAN_BUILD=false
CLEAN_MFEM=false
CLEAN_PLUGIN=false
UPDATE_IFSI=false
BUILD_TYPE="Release"

while [[ $# -gt 0 ]]; do
    case $1 in
        --febio-dir)
            FEBio_DIR="$2"
            shift 2
            ;;
        --build-type)
            BUILD_TYPE="$2"
            shift 2
            ;;
        --gslib-cflags)
            GSLIB_CFLAGS="$2"
            shift 2
            ;;
        --hypre-cflags)
            HYPRE_CFLAGS="$2"
            shift 2
            ;;
        --hypre-cxxflags)
            HYPRE_CXXFLAGS="$2"
            shift 2
            ;;
        --metis-cflags)
            METIS_CFLAGS="$2"
            shift 2
            ;;
        --metis-cxxflags)
            METIS_CXXFLAGS="$2"
            shift 2
            ;;
        --mfem-cxxflags)
            MFEM_CXXFLAGS="$2"
            shift 2
            ;;
        --cc)
            OVERRIDE_CC="$2"
            shift 2
            ;;
        --cxx)
            OVERRIDE_CXX="$2"
            shift 2
            ;;
        --apple-cc)
            APPLE_CC="$2"
            shift 2
            ;;
        --apple-cxx)
            APPLE_CXX="$2"
            shift 2
            ;;
        --mfem-options)
            MFEM_OPTIONS="$2"
            shift 2
            ;;
        --update-ifsi)
            UPDATE_IFSI=true
            shift
            ;;
        --clean)
            CLEAN_BUILD=true
            shift
            ;;
        --clean-mfem)
            CLEAN_MFEM=true
            shift
            ;;
        --clean-plugin)
            CLEAN_PLUGIN=true
            shift
            ;;
        --help)
            usage
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            usage
            ;;
    esac
done

# Validate required arguments
if [ -z "${FEBio_DIR}" ]; then
    echo -e "${RED}Error: --febio-dir is required${NC}"
    usage
fi

if [ ! -d "${FEBio_DIR}" ]; then
    echo -e "${RED}Error: FEBio installation directory not found at ${FEBio_DIR}${NC}"
    exit 1
fi

echo -e "${GREEN}=== FEBioMFEMiFSI Build Script ===${NC}"
echo "Project Root: ${PROJECT_ROOT}"
echo "Software Directory: ${SOFTWARE_DIR}"
echo "FEBio Directory: ${FEBio_DIR}"
echo "MFEM Directory: ${MFEM_DIR}"
echo "MFEM Build Directory: ${MFEM_BUILD_DIR}"
echo "Build Type: ${BUILD_TYPE}"

# Set up compiler args for CMake after parsing options
if [[ "${OSTYPE}" == "darwin"* ]]; then
    ARCH=$(uname -m)
    if [ "${ARCH}" != "x86_64" ]; then
        echo -e "${RED}Unsupported Apple architecture for this Plugin: ${ARCH}${NC}"
        exit 1
    fi
    # Prefer general overrides; fall back to apple-specific; default to clang/clang++
    CC_CHOICE="${OVERRIDE_CC}"
    CXX_CHOICE="${OVERRIDE_CXX}"
    if [ -z "${CC_CHOICE}" ]; then CC_CHOICE="${APPLE_CC}"; fi
    if [ -z "${CXX_CHOICE}" ]; then CXX_CHOICE="${APPLE_CXX}"; fi
    if [ -z "${CC_CHOICE}" ]; then CC_CHOICE="clang"; fi
    if [ -z "${CXX_CHOICE}" ]; then CXX_CHOICE="clang++"; fi
    # Validate chosen compilers
    if ! command -v "${CC_CHOICE}" >/dev/null 2>&1; then
        echo -e "${RED}Specified Apple C compiler not found: ${CC_CHOICE}${NC}"
        exit 1
    fi
    if ! command -v "${CXX_CHOICE}" >/dev/null 2>&1; then
        echo -e "${RED}Specified Apple C++ compiler not found: ${CXX_CHOICE}${NC}"
        exit 1
    fi
    echo -e "${YELLOW}Apple x86_64; using ${CC_CHOICE}/${CXX_CHOICE} for CMake builds${NC}"
    CMAKE_COMPILER_ARGS="-DCMAKE_C_COMPILER=${CC_CHOICE} -DCMAKE_CXX_COMPILER=${CXX_CHOICE}"
    MFEM_COMPILER_ARGS="-DCMAKE_CXX_COMPILER=${CXX_CHOICE}"
else
    # Non-Apple: apply overrides if provided
    # Record choices for use by make/autotools builds
    CC_CHOICE="${OVERRIDE_CC}"
    CXX_CHOICE="${OVERRIDE_CXX}"
    if [ -n "${OVERRIDE_CC}" ]; then
        if ! command -v "${OVERRIDE_CC}" >/dev/null 2>&1; then
            echo -e "${RED}Specified C compiler not found: ${OVERRIDE_CC}${NC}"
            exit 1
        fi
        CMAKE_COMPILER_ARGS="${CMAKE_COMPILER_ARGS} -DCMAKE_C_COMPILER=${OVERRIDE_CC}"
    fi
    if [ -n "${OVERRIDE_CXX}" ]; then
        if ! command -v "${OVERRIDE_CXX}" >/dev/null 2>&1; then
            echo -e "${RED}Specified C++ compiler not found: ${OVERRIDE_CXX}${NC}"
            exit 1
        fi
        CMAKE_COMPILER_ARGS="${CMAKE_COMPILER_ARGS} -DCMAKE_CXX_COMPILER=${OVERRIDE_CXX}"
        MFEM_COMPILER_ARGS="${MFEM_COMPILER_ARGS} -DCMAKE_CXX_COMPILER=${OVERRIDE_CXX}"
    fi
    if [ -n "${CMAKE_COMPILER_ARGS}" ]; then
        echo -e "${YELLOW}Using custom compilers for CMake builds: ${OVERRIDE_CC}/${OVERRIDE_CXX}${NC}"
    fi
fi

# Fetch MFEM if not present
# If building from scratch, we need to build dependencies first
if [ ! -d "${MFEM_DIR}" ]; then
    echo -e "${YELLOW}Building MFEM from scratch. Setting up dependencies...${NC}"
    
    # Create Software directory
    mkdir -p "${SOFTWARE_DIR}"
    
    # Build GSLIB
    if [ ! -d "${GSLIB_DIR}" ]; then
        echo -e "${YELLOW}Fetching GSLIB from GitHub...${NC}"
        git clone https://github.com/Nek5000/gslib.git "${GSLIB_DIR}"
        echo -e "${GREEN}GSLIB fetched successfully${NC}"
    else
        echo -e "${GREEN}GSLIB source already present${NC}"
    fi
    
    if [ ! -d "${GSLIB_BUILD_DIR}" ]; then
        echo -e "${YELLOW}Building GSLIB...${NC}"
        cd "${GSLIB_DIR}"
        if ! command -v mpicc >/dev/null 2>&1; then
            echo -e "${RED}Error: 'mpicc' not found. GSLIB requires MPI C compiler (mpicc). Please install MPI (e.g., OpenMPI) and ensure 'mpicc' is in PATH.${NC}"
            exit 1
        fi
        make CC="mpicc" CFLAGS="-fPIC ${GSLIB_CFLAGS}" -j${NUM_PROCS}
        mkdir -p "${GSLIB_BUILD_DIR}"
        echo -e "${GREEN}GSLIB built successfully${NC}"
    else
        echo -e "${GREEN}GSLIB already built${NC}"
    fi
    
    # Build HYPRE
    if [ ! -d "${HYPRE_DIR}" ]; then
        echo -e "${YELLOW}Fetching HYPRE from GitHub...${NC}"
        git clone https://github.com/hypre-space/hypre.git "${HYPRE_DIR}"
        echo -e "${GREEN}HYPRE fetched successfully${NC}"
    else
        echo -e "${GREEN}HYPRE source already present${NC}"
    fi
    
        if [ ! -f "${HYPRE_DIR}/src/hypre/lib/libHYPRE.a" ]; then
        echo -e "${YELLOW}Building HYPRE...${NC}"
        cd "${HYPRE_DIR}/src"
        ./configure --disable-fortran --with-extra-CXXFLAGS="-fPIC ${HYPRE_CXXFLAGS}" --with-extra-CFLAGS="-fPIC ${HYPRE_CFLAGS}"
        make -j${NUM_PROCS}
        echo -e "${GREEN}HYPRE built successfully${NC}"
    else
        echo -e "${GREEN}HYPRE already built${NC}"
    fi
    
    # Build METIS (v5.1.0)
    if [ ! -L "${METIS_DIR}" ]; then
        echo -e "${YELLOW}Fetching METIS from MFEM tpls...${NC}"
        
        # METIS v5.1.0 archive and extracted directory name
        METIS_ARCHIVE="metis-5.1.0.tar.gz"
        METIS_EXTRACTED_DIR="metis-5.1.0"
        
        cd "${SOFTWARE_DIR}"
        echo -e "${YELLOW}Downloading ${METIS_ARCHIVE}...${NC}"
        
        # Download from MFEM tpls repository (use raw.githubusercontent to avoid HTML redirect)
        DOWNLOAD_URL="https://raw.githubusercontent.com/mfem/tpls/master/${METIS_ARCHIVE}"
        echo "  Download URL: ${DOWNLOAD_URL}"
        
        if ! curl -L -f -o "${METIS_ARCHIVE}" "${DOWNLOAD_URL}"; then
            echo -e "${YELLOW}Primary download failed, trying fallback URL...${NC}"
            DOWNLOAD_URL="https://github.com/mfem/tpls/raw/master/${METIS_ARCHIVE}"
            echo "  Fallback URL: ${DOWNLOAD_URL}"
            curl -L -f -o "${METIS_ARCHIVE}" "${DOWNLOAD_URL}" || {
                echo -e "${RED}Error: Failed to download METIS from MFEM tpls repository${NC}"
                exit 1
            }
        fi
        
        # Check if download was successful by verifying file integrity
        if ! tar -tzf "${METIS_ARCHIVE}" > /dev/null 2>&1; then
            echo -e "${RED}Error: Downloaded file is not a valid tar.gz archive${NC}"
            rm -f "${METIS_ARCHIVE}"
            exit 1
        fi
        
        tar -xzf "${METIS_ARCHIVE}"
        rm "${METIS_ARCHIVE}"
        
        # Create/update symlink from metis-X.X.X to metis
        if [ -d "${METIS_EXTRACTED_DIR}" ]; then
            ln -sfn "${METIS_EXTRACTED_DIR}" metis
            echo -e "${GREEN}METIS fetched and symlinked successfully${NC}"
        else
            echo -e "${RED}Error: METIS extraction failed - directory ${METIS_EXTRACTED_DIR} not found${NC}"
            exit 1
        fi
    else
        echo -e "${GREEN}METIS source already present${NC}"
    fi
    
    if [ ! -f "${METIS_DIR}/lib/libmetis.a" ]; then
        echo -e "${YELLOW}Building METIS 5.1.0...${NC}"
        cd "${METIS_DIR}"
        # Configure METIS with CMake directly to ensure flags are applied
        METIS_CMAKE_CMD="cmake ${CMAKE_COMPILER_ARGS} -S '${METIS_DIR}' -B '${METIS_DIR}/lib' -DCMAKE_VERBOSE_MAKEFILE=1 -DGKLIB_PATH='${METIS_DIR}/GKlib' -DCMAKE_C_FLAGS='-fPIC ${METIS_CFLAGS}' -DCMAKE_CXX_FLAGS='-fPIC ${METIS_CXXFLAGS}'"
        echo "Running: ${METIS_CMAKE_CMD}"
        eval "${METIS_CMAKE_CMD}"
        cmake --build "${METIS_DIR}/lib" -j${NUM_PROCS}
        # Copy archive to the lib directory root for MFEM to find at METIS_DIR/lib
        if [ -f "${METIS_DIR}/lib/libmetis/libmetis.a" ]; then
            cp "${METIS_DIR}/lib/libmetis/libmetis.a" "${METIS_DIR}/lib/"
        fi
        echo -e "${GREEN}METIS built successfully${NC}"
    else
        echo -e "${GREEN}METIS already built${NC}"
    fi
    
    # Fetch MFEM
    echo -e "${YELLOW}Fetching MFEM from GitHub...${NC}"
    git clone https://github.com/mfem/mfem.git "${MFEM_DIR}"
    echo -e "${GREEN}MFEM fetched successfully${NC}"
else
    echo -e "${GREEN}MFEM source already present at ${MFEM_DIR}${NC}"
fi

# Check if MFEM source exists
if [ ! -d "${MFEM_DIR}" ]; then
    echo -e "${RED}Error: MFEM directory not found at ${MFEM_DIR}${NC}"
    exit 1
fi

# Build MFEM if needed
if [ ! -d "${MFEM_BUILD_DIR}" ] || [ "$CLEAN_BUILD" = true ] || [ "$CLEAN_MFEM" = true ]; then
    echo -e "${YELLOW}Building MFEM...${NC}"
    
    if { [ "$CLEAN_BUILD" = true ] || [ "$CLEAN_MFEM" = true ]; } && [ -d "${MFEM_BUILD_DIR}" ]; then
        rm -rf "${MFEM_BUILD_DIR}"
    fi
    
    mkdir -p "${MFEM_BUILD_DIR}"
    cd "${MFEM_BUILD_DIR}"
    
    # CMake configuration command for MFEM
    cmake_cmd="cmake ${MFEM_COMPILER_ARGS} -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DCMAKE_CXX_FLAGS='-fPIC ${MFEM_CXXFLAGS}' -DMFEM_USE_GSLIB=ON -DGSLIB_DIR=${GSLIB_BUILD_DIR} -DMFEM_USE_MPI=ON -DHYPRE_DIR=${HYPRE_BUILD_DIR} -DMFEM_USE_METIS=ON -DMETIS_DIR=${METIS_BUILD_DIR} .. ${MFEM_OPTIONS}"
    echo "Running: ${cmake_cmd}"
    eval "${cmake_cmd}"
    
    make -j${NUM_PROCS}
    
    echo -e "${GREEN}MFEM built successfully${NC}"
else
    echo -e "${GREEN}MFEM already built${NC}"
fi

# Update git subtree (MFEMFSI)
cd "${PROJECT_ROOT}"

# Check if MFEMFSI directory exists; always add if missing
if [ ! -d "${PROJECT_ROOT}/MFEMFSI" ]; then
    if git status --short | grep -q .; then
        echo -e "${RED}Error: Cannot add MFEMFSI subtree with uncommitted changes. Please commit or stash your changes first.${NC}"
        exit 1
    fi
    
    echo -e "${YELLOW}Adding MFEMFSI subtree...${NC}"
    
    # Check if iFSI remote exists
    if ! git remote get-url iFSI > /dev/null 2>&1; then
        echo -e "${YELLOW}Adding iFSI remote...${NC}"
        git remote add iFSI https://github.com/ryanblack129/iFSI.git
        git fetch iFSI
        echo -e "${GREEN}iFSI remote added${NC}"
    fi
    
    git subtree add --prefix=MFEMFSI iFSI main --squash
    echo -e "${GREEN}MFEMFSI subtree added${NC}"
elif [ "$UPDATE_IFSI" = true ]; then
    # Only update if explicitly requested via --update-ifsi
    if git status --short | grep -q .; then
        echo -e "${RED}Error: Cannot update iFSI subtree with uncommitted changes. Please commit or stash your changes first.${NC}"
        exit 1
    fi
    
    echo -e "${YELLOW}Updating MFEMFSI subtree...${NC}"
    # Check if iFSI remote exists
    if ! git remote get-url iFSI > /dev/null 2>&1; then
        echo -e "${YELLOW}Adding iFSI remote...${NC}"
        git remote add iFSI https://github.com/ryanblack129/iFSI.git
        git fetch iFSI
        echo -e "${GREEN}iFSI remote added${NC}"
    fi
    
    git subtree pull --prefix=MFEMFSI iFSI main --squash
    echo -e "${GREEN}MFEMFSI subtree updated${NC}"
else
    echo -e "${GREEN}Skipping iFSI subtree update (use --update-ifsi to update)${NC}"
fi

# Build FEBioMFEMiFSI project
echo -e "${YELLOW}Building FEBioMFEMiFSI...${NC}"

# Create build directory if it doesn't exist
if [ ! -d "${PROJECT_ROOT}/build" ]; then
    mkdir -p "${PROJECT_ROOT}/build"
fi

if [ "$CLEAN_BUILD" = true ] || [ "$CLEAN_PLUGIN" = true ]; then
    rm -rf "${PROJECT_ROOT}/build"
    mkdir -p "${PROJECT_ROOT}/build"
fi

cd "${PROJECT_ROOT}/build"

# CMake configuration for FEBioMFEMiFSI
# Set MFEM_DIR to point to the MFEM build directory and FEBio_DIR
cmake_cmd="cmake ${CMAKE_COMPILER_ARGS} -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DMFEM_DIR=${MFEM_BUILD_DIR} -DFEBio_DIR=${FEBio_DIR} .."
echo "Running: ${cmake_cmd}"
eval "${cmake_cmd}"

# Build
make -j${NUM_PROCS}

echo -e "${GREEN}=== Build completed successfully ===${NC}"
echo "Library location: ${PROJECT_ROOT}/build/lib/"
