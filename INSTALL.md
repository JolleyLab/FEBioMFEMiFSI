# FEBioMFEMiFSI Build Instructions

This document provides instructions on how to build the FEBioMFEMiFSI project using the provided `build.sh` script.

## Overview

The `build.sh` script is a comprehensive build automation tool that:
1. Fetches and builds all required dependencies (GSLIB, HYPRE, METIS, MFEM)
2. Updates the git subtree for the iFSI code
3. Builds the FEBioMFEMiFSI plugin

## Prerequisites

Before running the build script, ensure you have:

- **CMake** (version 3.10 or higher)
- **C/C++ Compiler** (gcc/clang with C++17 support)
- **Make** or compatible build system
- **Git** (for fetching dependencies and managing subtrees)
- **FEBio SDK** (FEBio installation)
- **MPI** (for parallel computing support)

### macOS support
- Supported only on Intel macOS (x86_64). Apple Silicon (arm64) is not supported for this build flow.
- On macOS, the script enforces x86_64 and configures CMake-based builds with `clang/clang++` by default. You can override compilers using `--cc` and `--cxx` (e.g., `--cc gcc-13 --cxx g++-13`).

## Quick Start

The simplest way to build is:

```bash
./build.sh --febio-sdk /path/to/FEBio
```

## Build Script Options

### Required Options

- `--febio-sdk PATH` - Path to FEBio SDK root directory (REQUIRED)
  - Example: `--febio-sdk /Users/user/Software/FEBio`

### Optional Options

- `--build-type TYPE` - Build type: `Debug` or `Release` (default: `Release`)
- `--gslib-cflags FLAGS` - Additional C compiler flags for GSLIB
- `--hypre-cflags FLAGS` - Additional C compiler flags for HYPRE
- `--hypre-cxxflags FLAGS` - Additional C++ compiler flags for HYPRE
- `--metis-cflags FLAGS` - Additional C compiler flags for METIS
- `--metis-cxxflags FLAGS` - Additional C++ compiler flags for METIS
- `--mfem-cxxflags FLAGS` - Additional C++ compiler flags for MFEM
 - `--cc CC` - Override C compiler for CMake builds (macOS defaults to clang)
 - `--cxx CXX` - Override C++ compiler for CMake builds (macOS defaults to clang++)
- `--mfem-options STR` - Additional CMake options for MFEM
  - Useful for enabling/disabling specific MFEM features
  - Example: `--mfem-options '-DMFEM_USE_MPI=ON'`
- `--update-ifsi` - Update iFSI subtree from remote repository
  - Requires a clean working tree (no uncommitted changes)
  - Skipped by default; only runs when this flag is provided
  - The subtree is automatically added on first build if missing
- `--clean` - Clean and rebuild both MFEM and FEBioMFEMiFSI plugin
  - Does not clean dependencies (GSLIB, HYPRE, METIS)
- `--clean-mfem` - Clean and rebuild MFEM only
- `--clean-plugin` - Clean and rebuild FEBioMFEMiFSI plugin only
- `--help` - Display help message and exit

## Usage Examples

### Basic Release Build
```bash
./build.sh --febio-sdk /path/to/FEBio
```

### Clean Rebuild
```bash
# Clean and rebuild everything
./build.sh --febio-sdk /path/to/FEBio --clean

# Clean and rebuild MFEM only
./build.sh --febio-sdk /path/to/FEBio --clean-mfem

# Clean and rebuild plugin only
./build.sh --febio-sdk /path/to/FEBio --clean-plugin
```

### Update iFSI Subtree
```bash
./build.sh --febio-sdk /path/to/FEBio --update-ifsi
```

## Build Output

After successful build, you'll find:
- **Plugin library**: `./build/lib/libFEBioMFEMiFSI.so` (Linux/macOS)
- **Plugin library**: `./build/lib/FEBioMFEMiFSI.dylib` (macOS)

The script creates the following build directory structure:

```
.
├── Software/                      # Third-party libraries and dependencies
│   ├── mfem/                      # MFEM source and build
│   │   ├── build-mfem/            # MFEM build directory
│   │   └── [MFEM source files]
│   ├── gslib/                     # GSLIB source and build
│   │   ├── build/                 # GSLIB build directory
│   │   └── [GSLIB source files]
│   ├── hypre/                     # HYPRE source and build (autotools, in-place)
│   │   ├── src/                   # HYPRE source with libHYPRE.a built in-place
│   │   └── [HYPRE configuration]
│   ├── metis-5.1.0/               # METIS 5 extracted source
│   │   ├── include/               # METIS headers (metis.h)
│   │   ├── lib/                   # Built METIS static library location
│   │   │   └── libmetis.a
│   │   └── [METIS source files]
│   └── metis → metis-5.1.0        # Symlink to current METIS version
├── build/                         # FEBioMFEMiFSI build
│   ├── lib/                       # Output library
│   └── CMakeFiles/
├── MFEMFSI/                       # iFSI subtree (git managed)
│   └── [iFSI source files]
└── [Other project files]
```

---

## Manual Build Instructions

If the automated `build.sh` script fails or you need more control over the build process, follow these manual steps to build each dependency and the plugin.

### Environment Setup

Set up directory structure and environment variables:

```bash
export PROJECT_ROOT="$(pwd)"
export SOFTWARE_DIR="${PROJECT_ROOT}/Software"
export FEBio_SDK="/path/to/your/FEBio"
mkdir -p "${SOFTWARE_DIR}"
```

### 1. Build GSLIB

GSLIB is built using Make with the MPI C compiler.

**Default Build Options:**
- Compiler: `mpicc` (required)
- Flags: `-fPIC` (mandatory for shared library integration)

**Build Steps:**

```bash
cd "${SOFTWARE_DIR}"
git clone https://github.com/Nek5000/gslib.git
cd gslib
make CC=mpicc CFLAGS="-fPIC"
mkdir -p build
```

**Output:**
- Built library files in `${SOFTWARE_DIR}/gslib/`
- Build marker directory: `${SOFTWARE_DIR}/gslib/build`

### 2. Build HYPRE

HYPRE is built using autotools with an in-place build in the `src/` directory.

**Default Build Options:**
- Build system: autotools (configure + make)
- Fortran: disabled (`--disable-fortran`)
- Flags: `-fPIC` (mandatory)
- Compiler: System defaults (not overridden)

**Build Steps:**

```bash
cd "${SOFTWARE_DIR}"
git clone https://github.com/hypre-space/hypre.git
cd hypre/src
./configure --disable-fortran \
            --with-extra-CFLAGS="-fPIC" \
            --with-extra-CXXFLAGS="-fPIC"
make
```

**Output:**
- Static library: `${SOFTWARE_DIR}/hypre/src/hypre/lib/libHYPRE.a`
- MFEM will use: `HYPRE_DIR=${SOFTWARE_DIR}/hypre/src/hypre/lib`

### 3. Build METIS 5.1.0

METIS is built using CMake with manual download from the MFEM tpls repository.

**Default Build Options:**
- Version: 5.1.0 (downloaded from MFEM tpls)
- Build system: CMake
- Flags: `-fPIC` (mandatory)
- Compiler: CMake defaults or user-specified via `--cc`/`--cxx`
- Build directory: `lib/` subdirectory within METIS source

**Build Steps:**

```bash
cd "${SOFTWARE_DIR}"

# Download METIS 5.1.0
curl -L -f -o metis-5.1.0.tar.gz \
  https://raw.githubusercontent.com/mfem/tpls/master/metis-5.1.0.tar.gz

# Verify download integrity
tar -tzf metis-5.1.0.tar.gz > /dev/null || { echo "Invalid archive"; exit 1; }

# Extract and create symlink
tar -xzf metis-5.1.0.tar.gz
rm metis-5.1.0.tar.gz
ln -sfn metis-5.1.0 metis

# Build METIS
cd metis
cmake -S . -B lib \
      -DCMAKE_VERBOSE_MAKEFILE=1 \
      -DGKLIB_PATH="$(pwd)/GKlib" \
      -DCMAKE_C_FLAGS="-fPIC" \
      -DCMAKE_CXX_FLAGS="-fPIC"
cmake --build lib

# Copy library to expected location
cp lib/libmetis/libmetis.a lib/
```

**Output:**
- Static library: `${SOFTWARE_DIR}/metis/lib/libmetis.a`
- Headers: `${SOFTWARE_DIR}/metis/include/metis.h`
- MFEM will use: `METIS_DIR=${SOFTWARE_DIR}/metis`

### 4. Build MFEM

MFEM is built using CMake with integrations for GSLIB, HYPRE, and METIS.

**Default Build Options:**
- Build system: CMake
- Build type: Release (or Debug via `--build-type`)
- Flags: `-fPIC` (mandatory)
- Compiler: On macOS x86_64, defaults to `clang++`; otherwise CMake defaults or user-specified
- Integrations enabled: GSLIB, HYPRE, METIS, MPI
- Build directory: `build/` subdirectory

**Build Steps:**

```bash
cd "${SOFTWARE_DIR}"
git clone https://github.com/mfem/mfem.git
cd mfem
mkdir -p build
cd build

# Configure MFEM with dependencies
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-fPIC" \
  -DMFEM_USE_GSLIB=ON \
  -DGSLIB_DIR="${SOFTWARE_DIR}/gslib/build" \
  -DMFEM_USE_MPI=ON \
  -DHYPRE_DIR="${SOFTWARE_DIR}/hypre/src/hypre/lib" \
  -DMFEM_USE_METIS=ON \
  -DMETIS_DIR="${SOFTWARE_DIR}/metis"

# Build
make
```

**macOS-specific compiler override:**

If you need to override compilers on macOS (or any platform):

```bash
cmake .. \
  -DCMAKE_CXX_COMPILER=clang++ \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-fPIC" \
  -DMFEM_USE_GSLIB=ON \
  -DGSLIB_DIR="${SOFTWARE_DIR}/gslib/build" \
  -DMFEM_USE_MPI=ON \
  -DHYPRE_DIR="${SOFTWARE_DIR}/hypre/src/hypre/lib" \
  -DMFEM_USE_METIS=ON \
  -DMETIS_DIR="${SOFTWARE_DIR}/metis"
```

**Output:**
- Build directory: `${SOFTWARE_DIR}/mfem/build-mfem/`
- MFEM config: `${SOFTWARE_DIR}/mfem/build-mfem/MFEMConfig.cmake`

### 5. Update iFSI Subtree (Optional)

If you need to update the MFEMFSI subtree and have a clean working tree:

```bash
cd "${PROJECT_ROOT}"

# Add remote if not already added
git remote get-url iFSI || git remote add iFSI https://github.com/ryanblack129/iFSI.git
git fetch iFSI

# Add or update subtree
if [ ! -d MFEMFSI ]; then
  git subtree add --prefix=MFEMFSI iFSI main --squash
else
  git subtree pull --prefix=MFEMFSI iFSI main --squash
fi
```

**Note:** Git subtree operations require a clean working tree (no uncommitted changes).

### 6. Build FEBioMFEMiFSI Plugin

Build the final plugin linking against FEBio SDK and MFEM with all dependencies.

**Default Build Options:**
- Build system: CMake
- Build type: Release (or Debug)
- Compiler: Matches MFEM compiler choice (CMake defaults or user-specified)
- Build directory: `build/`

**Build Steps:**

```bash
cd "${PROJECT_ROOT}"
mkdir -p build
cd build

# Configure plugin
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DMFEM_DIR="${SOFTWARE_DIR}/mfem/build-mfem" \
  -DFEBio_SDK="${FEBio_SDK}"

# Build
make
```

**macOS-specific compiler override:**

```bash
cmake .. \
  -DCMAKE_C_COMPILER=clang \
  -DCMAKE_CXX_COMPILER=clang++ \
  -DCMAKE_BUILD_TYPE=Release \
  -DMFEM_DIR="${SOFTWARE_DIR}/mfem/build-mfem" \
  -DFEBio_SDK="${FEBio_SDK}"
```

**Output:**
- Plugin library: `${PROJECT_ROOT}/build/lib/libFEBioMFEMiFSI.so` (Linux)
- Plugin library: `${PROJECT_ROOT}/build/lib/FEBioMFEMiFSI.dylib` (macOS)

### Clean Rebuild

To perform a clean rebuild of MFEM and the plugin (without rebuilding dependencies):

```bash
rm -rf "${SOFTWARE_DIR}/mfem/build-mfem"
rm -rf "${PROJECT_ROOT}/build"
# Then repeat steps 4 and 6
```

To rebuild everything from scratch:

```bash
rm -rf "${SOFTWARE_DIR}"
rm -rf "${PROJECT_ROOT}/build"
# Then repeat all steps 1-6
```
