# Building and Installing MDcraft

This README focuses on practical build and installation workflows for MDcraft on Linux. The recommended entry point is the automated Linux bootstrap script, but the traditional `cmake` and `ccmake` workflow is still fully supported and documented below.

There are two supported Linux workflows:

- `scripts/bootstrap_linux.sh` for the common Python-enabled workflow
- direct `cmake` or `ccmake` configuration for manual, offline, or highly customized builds

If you are starting from scratch, read the Linux quick start first. If you prefer to work directly with `cmake` or `ccmake`, jump to [Linux manual start](#linux-manual-start).

## Linux quick start

### Recommended workflow: `scripts/bootstrap_linux.sh`

Create and activate the Python environment where you want to use `mdcraft`, then run:

```bash
./scripts/bootstrap_linux.sh
```

This script is the main entry point for Linux users. It is usually easier than configuring the project manually because it installs Python requirements, configures CMake, builds the project, and installs it in one pass.

Current defaults:

- CMake executable: `/usr/bin/cmake`
- C and C++ compilers: `/usr/bin/cc` and `/usr/bin/c++`
- Python interpreter: the current `python3`, unless `VENV_DIR` or `PYTHON_BIN` is set
- build directory: `build/linux-release`
- install prefix: `.local`
- Python install layout: `site-packages`
- MPI: enabled by default
- Python bindings: enabled
- tests: enabled
- example problems: enabled
- automatic download of missing CMake dependencies: enabled
- thread implementation: `STD`

By default the script installs Python packages from [python/requirements.txt](python/requirements.txt), which currently includes:

- `numpy`
- `scipy`
- `matplotlib`
- `h5py`
- `pybind11`
- `ase`

When MPI is enabled, bootstrap also installs `mpi4py` for the selected Python interpreter. If a dependency required by the selected options is missing, CMake stops with a clear error during configuration.

After installation with the default `site-packages` layout, you should be able to run Python scripts directly from the selected environment without setting `PYTHONPATH` or `LD_LIBRARY_PATH` manually.

### `bootstrap_linux.sh --help`

The script has a built-in help page:

```bash
./scripts/bootstrap_linux.sh --help
```

It prints the supported environment variables and current defaults. The most important variables are:

Toolchain and build layout:

- `CMAKE_BIN`: path to `cmake`
- `CMAKE_GENERATOR`: CMake generator such as `Unix Makefiles` or `Ninja`
- `BUILD_DIR`: CMake build directory
- `INSTALL_PREFIX`: value for `CMAKE_INSTALL_PREFIX`
- `BUILD_TYPE`: value for `CMAKE_BUILD_TYPE`
- `BUILD_PARALLELISM`: optional value for `cmake --build --parallel`
- `C_COMPILER`: value for `CMAKE_C_COMPILER`
- `CXX_COMPILER`: value for `CMAKE_CXX_COMPILER`

Python environment:

- `PYTHON_BIN`: Python interpreter used by bootstrap and passed to CMake as `Python_EXECUTABLE`
- `VENV_DIR`: create or reuse this virtual environment and use its `bin/python`
- `VENV_BASE_PYTHON`: interpreter used to create `VENV_DIR`
- `SKIP_PIP_INSTALL`: skip `pip` installation steps when set to `ON`
- `PYTHON_INSTALL_LAYOUT`: choose the Python installation layout used by CMake

Feature selection:

- `ENABLE_MPI`: maps to `mdcraft_ENABLE_MPI`
- `ENABLE_TESTS`: maps to `mdcraft_ENABLE_TESTS`
- `ENABLE_PROBLEMS`: maps to `mdcraft_ENABLE_PROBLEMS`
- `ENABLE_BPNN`: maps to `mdcraft_ENABLE_BPNN`
- `ENABLE_DEEPMD`: maps to `mdcraft_ENABLE_DeePMD`
- `ENABLE_MLIP4`: maps to `mdcraft_ENABLE_MLIP4`
- `THREADS_TYPE`: `STD` or `TBB`
- `FETCH_MISSING_DEPS`: maps to `mdcraft_FETCH_MISSING_DEPS`
- `CONFIGURE_ONLY`: stop after CMake configure

MPI-specific settings:

- `MPI_C_COMPILER`: passed to CMake as `MPI_C_COMPILER`
- `MPI_CXX_COMPILER`: passed to CMake as `MPI_CXX_COMPILER`
- `MPIEXEC_BIN`: passed to CMake as `MPIEXEC_EXECUTABLE`
- `MPI4PY_PIP_SPEC`: package spec used when bootstrap installs `mpi4py`

Additional command-line arguments after the script name are forwarded directly to the CMake configure command. This is the right place for explicit `-D...` options such as `-Dmdcraft_FETCH_MISSING_DEPS=OFF`.

### Common bootstrap examples

Default release build:

```bash
VENV_DIR=$PWD/.venv-mpi ./scripts/bootstrap_linux.sh
```

Non-MPI build:

```bash
VENV_DIR=$PWD/.venv-serial ENABLE_MPI=OFF ./scripts/bootstrap_linux.sh
```

Configure only, without downloads, using explicitly installed dependencies:

```bash
SKIP_PIP_INSTALL=ON \
CONFIGURE_ONLY=ON \
FETCH_MISSING_DEPS=OFF \
./scripts/bootstrap_linux.sh \
  -DEigen3_DIR=/opt/eigen/share/eigen3/cmake \
  -Dpybind11_DIR=/opt/pybind11/share/cmake/pybind11 \
  -DGTest_DIR=/opt/googletest/lib/cmake/GTest \
  -DTBB_DIR=/opt/tbb/lib/cmake/TBB
```

When rebuilding, it is safer to start from a fresh build directory, or remove `CMakeCache.txt` and `CMakeFiles/`, when changing:

- the CMake generator
- the C and C++ compilers
- the Python interpreter or virtual environment
- the MPI implementation

If a previous configure step downloaded dependencies into `build/.../_deps/` and you now want a strictly local or offline build, a clean build directory is also recommended.

### MPI with a local virtual environment

For an isolated local environment:

```bash
VENV_DIR=$PWD/.venv-mpi ENABLE_MPI=ON ./scripts/bootstrap_linux.sh
```

In this mode the script:

- creates the virtual environment if it does not exist
- installs Python requirements into that environment
- installs `mpi4py` using the selected MPI compiler wrapper
- passes that exact interpreter to CMake as `Python_EXECUTABLE`
- passes the selected MPI wrappers to CMake

After that you can run the MPI example test directly from the virtual environment:

```bash
source ./.venv-mpi/bin/activate
mpirun -np 2 python test/python/solver/run_ball_ar.py
```

### Manual CMake configuration

You can skip bootstrap and configure MDcraft directly with CMake. The project supports installing Python modules into the selected interpreter's `site-packages` directory.

Example with MPI:

```bash
python3 -m venv .venv-mpi
source .venv-mpi/bin/activate
pip install -r python/requirements.txt
pip install mpi4py

cmake -S . -B build/linux-mpi-release -G "Unix Makefiles" \
  -DCMAKE_C_COMPILER=/usr/bin/cc \
  -DCMAKE_CXX_COMPILER=/usr/bin/c++ \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$PWD/.local \
  -DPython_EXECUTABLE=$PWD/.venv-mpi/bin/python \
  -DMPI_C_COMPILER=/usr/bin/mpicc \
  -DMPI_CXX_COMPILER=/usr/bin/mpicxx \
  -DMPIEXEC_EXECUTABLE=/usr/bin/mpirun \
  -Dmdcraft_ENABLE_MPI=ON \
  -Dmdcraft_ENABLE_PYTHON=ON \
  -Dmdcraft_ENABLE_TESTS=ON \
  -Dmdcraft_ENABLE_PROBLEMS=ON \
  -Dmdcraft_THREADS_TYPE=STD \
  -Dmdcraft_FETCH_MISSING_DEPS=ON \
  -Dmdcraft_PYTHON_INSTALL_LAYOUT=site-packages

cmake --build build/linux-mpi-release --parallel
cmake --install build/linux-mpi-release
```

### Reconfiguring with `ccmake`

You can change cache values interactively:

```bash
cd build/linux-release
ccmake .
```

or from the repository root:

```bash
ccmake build/linux-release
```

Inside `ccmake`:

1. press `t` to show advanced cache variables when needed
2. edit `mdcraft_ENABLE_*`, `mdcraft_THREADS_TYPE`, `mdcraft_FETCH_MISSING_DEPS`, `Python_EXECUTABLE`, `Eigen3_DIR`, `pybind11_DIR`, `GTest_DIR`, `TBB_DIR`, `DeePMD_DIR`, `mlip4_DIR` and other relevant entries
3. press `c` to reconfigure
4. press `g` to regenerate the build system

Then rebuild and reinstall:

```bash
cmake --build build/linux-release --parallel
cmake --install build/linux-release
```

## Linux manual start

This section summarizes the direct `cmake` and `ccmake` workflow without the bootstrap script.

### Prerequisites

For a typical Linux build you need:

- `cmake` and `ccmake`
- a C and C++ compiler
- Python 3 with development headers
- Python packages from `python/requirements.txt`
- `mpi4py` if `mdcraft_ENABLE_MPI=ON`
- a local installation of any optional package you enable, such as DeePMD

Examples on Debian or Ubuntu:

```bash
sudo apt install cmake-curses-gui python3 python3-dev python3-venv python3-pip libopenmpi-dev
python3 -m venv .mdenv
source .mdenv/bin/activate
python -m pip install -U pip
python -m pip install -r python/requirements.txt
python -m pip install mpi4py
```

### Common package locations

If you install dependencies locally, CMake usually needs the directory that contains the package config file:

- `Eigen3_DIR` -> directory containing `Eigen3Config.cmake`
- `pybind11_DIR` -> directory containing `pybind11Config.cmake`
- `GTest_DIR` -> directory containing `GTestConfig.cmake`
- `TBB_DIR` -> directory containing `TBBConfig.cmake`
- `DeePMD_DIR` -> directory containing `DeePMDConfig.cmake`

Examples:

```text
/opt/eigen/share/eigen3/cmake
/opt/pybind11/share/cmake/pybind11
/opt/googletest/lib/cmake/GTest
/opt/tbb/lib/cmake/TBB
/home/username/soft/install/deepmd/lib/cmake/DeePMD
```

### Direct `cmake` workflow

MPI-enabled example:

```bash
cmake -S . -B build/linux-mpi-release -G "Unix Makefiles" \
  -DCMAKE_C_COMPILER=/usr/bin/cc \
  -DCMAKE_CXX_COMPILER=/usr/bin/c++ \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$PWD/.local \
  -DPython_EXECUTABLE=$PWD/.mdenv/bin/python \
  -DMPI_C_COMPILER=/usr/bin/mpicc \
  -DMPI_CXX_COMPILER=/usr/bin/mpicxx \
  -DMPIEXEC_EXECUTABLE=/usr/bin/mpirun \
  -Dmdcraft_ENABLE_MPI=ON \
  -Dmdcraft_ENABLE_PYTHON=ON \
  -Dmdcraft_ENABLE_TESTS=ON \
  -Dmdcraft_ENABLE_PROBLEMS=ON

cmake --build build/linux-mpi-release --parallel
cmake --install build/linux-mpi-release
```

If you need local-only or offline configuration, add `-Dmdcraft_FETCH_MISSING_DEPS=OFF` and the relevant `*_DIR` variables.

### Interactive `ccmake` workflow

```bash
mkdir mdcraft-build
cd mdcraft-build
ccmake /path/to/src/mdcraft
```

Key cache entries to review:

```text
CMAKE_BUILD_TYPE
CMAKE_C_COMPILER
CMAKE_CXX_COMPILER
CMAKE_INSTALL_PREFIX
Python_EXECUTABLE
Eigen3_DIR
pybind11_DIR
GTest_DIR
TBB_DIR
DeePMD_DIR
mlip4_DIR
mdcraft_ENABLE_MPI
mdcraft_ENABLE_PYTHON
mdcraft_ENABLE_TESTS
mdcraft_ENABLE_PROBLEMS
mdcraft_ENABLE_BPNN
mdcraft_ENABLE_DeePMD
mdcraft_ENABLE_MLIP4
mdcraft_FETCH_MISSING_DEPS
mdcraft_THREADS_TYPE
```

Inside `ccmake`:

1. press `c` for the initial configure step
2. press `t` if you want to see advanced cache variables
3. edit the relevant options and dependency paths
4. press `c` again to reconfigure
5. press `g` to generate the build system

Then run:

```bash
cmake --build . --parallel
cmake --install .
ctest --output-on-failure
```

### Rebuild after changing options

You can usually keep the same build directory when changing feature flags or dependency paths. It is safer to recreate the build directory when changing the generator, compiler, Python interpreter, or MPI implementation.

### Running after installation

With the default `site-packages` layout, activate the target Python environment and run your scripts directly:

```bash
source path/to/env/.mdenv/bin/activate
mpirun -np 2 python test/python/solver/run_ball_ar.py
```

If you are using a custom manual layout, set `LD_LIBRARY_PATH` and `PYTHONPATH` as needed for that installation.

For example, if shared libraries are installed into `/home/username/soft/mdcraft/lib` and Python modules into `/home/username/soft/mdcraft/python`:

```bash
export LD_LIBRARY_PATH=/home/username/soft/mdcraft/lib:${LD_LIBRARY_PATH}
export PYTHONPATH=/home/username/soft/mdcraft/python:${PYTHONPATH}
```

On macOS use `DYLD_LIBRARY_PATH` instead of `LD_LIBRARY_PATH`:

```bash
export DYLD_LIBRARY_PATH=/Users/username/soft/mdcraft/lib:${DYLD_LIBRARY_PATH}
export PYTHONPATH=/Users/username/soft/mdcraft/python:${PYTHONPATH}
```
