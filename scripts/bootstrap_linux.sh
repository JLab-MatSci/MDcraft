#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CMAKE_BIN="${CMAKE_BIN:-/usr/bin/cmake}"
CMAKE_GENERATOR="${CMAKE_GENERATOR:-}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
VENV_DIR="${VENV_DIR:-}"
VENV_BASE_PYTHON="${VENV_BASE_PYTHON:-/usr/bin/python3}"
BUILD_DIR="${BUILD_DIR:-${ROOT_DIR}/build/linux-release}"
INSTALL_PREFIX="${INSTALL_PREFIX:-${ROOT_DIR}/.local}"
BUILD_TYPE="${BUILD_TYPE:-Release}"
BUILD_PARALLELISM="${BUILD_PARALLELISM:-}"
C_COMPILER="${C_COMPILER:-/usr/bin/cc}"
CXX_COMPILER="${CXX_COMPILER:-/usr/bin/c++}"
MPI_C_COMPILER="${MPI_C_COMPILER:-/usr/bin/mpicc}"
MPI_CXX_COMPILER="${MPI_CXX_COMPILER:-/usr/bin/mpicxx}"
MPIEXEC_BIN="${MPIEXEC_BIN:-/usr/bin/mpirun}"
MPI4PY_PIP_SPEC="${MPI4PY_PIP_SPEC:-mpi4py}"
ENABLE_MPI="${ENABLE_MPI:-ON}"
ENABLE_TESTS="${ENABLE_TESTS:-ON}"
ENABLE_PROBLEMS="${ENABLE_PROBLEMS:-ON}"
THREADS_TYPE="${THREADS_TYPE:-STD}"
ENABLE_BPNN="${ENABLE_BPNN:-OFF}"
FETCH_MISSING_DEPS="${FETCH_MISSING_DEPS:-ON}"
PYTHON_INSTALL_LAYOUT="${PYTHON_INSTALL_LAYOUT:-site-packages}"
SKIP_PIP_INSTALL="${SKIP_PIP_INSTALL:-OFF}"
if [[ -n "${ENABLE_DEEPMD+x}" && -n "${ENABLE_DeePMD+x}" && "${ENABLE_DEEPMD}" != "${ENABLE_DeePMD}" ]]; then
  echo "ENABLE_DEEPMD and ENABLE_DeePMD are both set but differ. Please use one value." >&2
  exit 1
fi
ENABLE_DEEPMD="${ENABLE_DEEPMD:-${ENABLE_DeePMD:-OFF}}"
ENABLE_MLIP4="${ENABLE_MLIP4:-OFF}"
CONFIGURE_ONLY="${CONFIGURE_ONLY:-OFF}"
GENERATOR_ARGS=()

show_help() {
  cat <<EOF
Usage:
  ./scripts/bootstrap_linux.sh [CMake configure args...]
  ./scripts/bootstrap_linux.sh --help

This script optionally creates a Python virtual environment, installs Python
requirements, configures MDcraft, then builds and installs it.

Environment variables:
  CMAKE_BIN             CMake executable (default: ${CMAKE_BIN})
  CMAKE_GENERATOR       CMake generator, e.g. "Unix Makefiles" or "Ninja"
  BUILD_DIR             CMake build directory (default: ${BUILD_DIR})
  INSTALL_PREFIX        CMake install prefix (default: ${INSTALL_PREFIX})
  BUILD_TYPE            CMAKE_BUILD_TYPE (default: ${BUILD_TYPE})
  BUILD_PARALLELISM     Value for 'cmake --build --parallel N' (default: CMake default)
  C_COMPILER            C compiler (default: ${C_COMPILER})
  CXX_COMPILER          C++ compiler (default: ${CXX_COMPILER})

  PYTHON_BIN            Python interpreter to use (default: ${PYTHON_BIN})
  VENV_DIR              Create/use this virtual environment instead of PYTHON_BIN
  VENV_BASE_PYTHON      Python used to create VENV_DIR (default: ${VENV_BASE_PYTHON})
  SKIP_PIP_INSTALL      ON/OFF, skip pip upgrade/install steps (default: ${SKIP_PIP_INSTALL})
  PYTHON_INSTALL_LAYOUT 'site-packages' or 'legacy' (default: ${PYTHON_INSTALL_LAYOUT})

  ENABLE_MPI            ON/OFF, enable MPI build (default: ${ENABLE_MPI})
  MPI_C_COMPILER        MPI C wrapper (default: ${MPI_C_COMPILER})
  MPI_CXX_COMPILER      MPI C++ wrapper (default: ${MPI_CXX_COMPILER})
  MPIEXEC_BIN           mpirun/mpiexec executable (default: ${MPIEXEC_BIN})
  MPI4PY_PIP_SPEC       Package spec used for mpi4py install (default: ${MPI4PY_PIP_SPEC})

  ENABLE_TESTS          ON/OFF, build tests (default: ${ENABLE_TESTS})
  ENABLE_PROBLEMS       ON/OFF, build example problems (default: ${ENABLE_PROBLEMS})
  ENABLE_BPNN           ON/OFF, enable BPNN support (default: ${ENABLE_BPNN})
  ENABLE_DEEPMD         ON/OFF, enable DeePMD support (default: ${ENABLE_DEEPMD})
  ENABLE_MLIP4          ON/OFF, enable MLIP4 support (default: ${ENABLE_MLIP4})
  THREADS_TYPE          STD or TBB (default: ${THREADS_TYPE})
  FETCH_MISSING_DEPS    ON/OFF for mdcraft_FETCH_MISSING_DEPS (default: ${FETCH_MISSING_DEPS})
  CONFIGURE_ONLY        ON/OFF, skip build/install after configure (default: ${CONFIGURE_ONLY})

Notes:
  - Additional command-line arguments are forwarded to the CMake configure step.
  - Use forwarded -D options to set explicit dependency paths or to disable downloads.

Examples:
  VENV_DIR=\$PWD/.venv-mpi ENABLE_MPI=ON ./scripts/bootstrap_linux.sh

  VENV_DIR=\$PWD/.venv-serial ENABLE_MPI=OFF ./scripts/bootstrap_linux.sh

  VENV_DIR=\$PWD/.venv-mpi SKIP_PIP_INSTALL=ON CONFIGURE_ONLY=ON FETCH_MISSING_DEPS=OFF \\
    ./scripts/bootstrap_linux.sh \\
    -DEigen3_DIR=/opt/eigen/share/eigen3/cmake \\
    -Dpybind11_DIR=/opt/pybind11/share/cmake/pybind11 \\
    -DGTest_DIR=/opt/googletest/lib/cmake/GTest
EOF
}

if [[ $# -gt 0 ]]; then
  case "${1}" in
    -h|--help|help)
      show_help
      exit 0
      ;;
  esac
fi

EXTRA_CMAKE_ARGS=("$@")

if [[ -n "${CMAKE_GENERATOR}" ]]; then
  GENERATOR_ARGS=(-G "${CMAKE_GENERATOR}")
fi

if [[ ! -x "${CMAKE_BIN}" ]]; then
  echo "CMake not found at ${CMAKE_BIN}" >&2
  exit 1
fi

if [[ -n "${VENV_DIR}" ]]; then
  if [[ ! -x "${VENV_DIR}/bin/python" ]]; then
    echo "Creating virtual environment in ${VENV_DIR}"
    "${VENV_BASE_PYTHON}" -m venv "${VENV_DIR}"
  fi
  PYTHON_BIN="${VENV_DIR}/bin/python"
fi

if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
  echo "Python interpreter '${PYTHON_BIN}' was not found in PATH." >&2
  exit 1
fi

RUN_ENV=(
  env
  CC=
  CXX=
  CPP=
  FC=
  F77=
  F90=
  CFLAGS=
  CXXFLAGS=
  CPPFLAGS=
  LDFLAGS=
  DEBUG_CFLAGS=
  DEBUG_CXXFLAGS=
  DEBUG_CPPFLAGS=
  CMAKE_PREFIX_PATH=
  PKG_CONFIG_PATH=
)
if [[ "${ENABLE_MPI}" == "ON" ]]; then
  RUN_ENV+=(
    PYTHONPATH=
    LD_LIBRARY_PATH=
    CONDA_PREFIX=
    CONDA_DEFAULT_ENV=
    CONDA_PROMPT_MODIFIER=
    CONDA_SHLVL=0
    CONDA_PYTHON_EXE=
    CONDA_BUILD_SYSROOT=
    _CONDA_PYTHON_SYSCONFIGDATA_NAME=
  )
fi

PYTHON_PACKAGES=(-r "${ROOT_DIR}/python/requirements.txt")
if [[ "${ENABLE_BPNN}" == "ON" ]]; then
  PYTHON_PACKAGES+=(ase)
fi

echo "Using Python: $("${PYTHON_BIN}" -c 'import sys; print(sys.executable)')"
if [[ "${SKIP_PIP_INSTALL}" == "ON" ]]; then
  echo "Skipping pip installation because SKIP_PIP_INSTALL=ON"
else
  echo "Installing Python dependencies into the active environment"
  "${PYTHON_BIN}" -m pip install --upgrade pip
  "${PYTHON_BIN}" -m pip install "${PYTHON_PACKAGES[@]}"
fi

if [[ "${ENABLE_MPI}" == "ON" ]]; then
  if [[ "${SKIP_PIP_INSTALL}" == "ON" ]]; then
    echo "Skipping mpi4py installation because SKIP_PIP_INSTALL=ON"
  else
    echo "Installing ${MPI4PY_PIP_SPEC} with MPICC=${MPI_C_COMPILER}"
    MPICC="${MPI_C_COMPILER}" "${PYTHON_BIN}" -m pip install "${MPI4PY_PIP_SPEC}"
  fi
fi

echo "Configuring MDcraft in ${BUILD_DIR}"
CONFIGURE_ARGS=(
  -S "${ROOT_DIR}"
  -B "${BUILD_DIR}"
  -DCMAKE_C_COMPILER="${C_COMPILER}"
  -DCMAKE_CXX_COMPILER="${CXX_COMPILER}"
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"
  -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}"
  -DPython_EXECUTABLE="${PYTHON_BIN}"
  -Dmdcraft_ENABLE_PYTHON=ON
  -Dmdcraft_ENABLE_MPI="${ENABLE_MPI}"
  -Dmdcraft_ENABLE_TESTS="${ENABLE_TESTS}"
  -Dmdcraft_ENABLE_PROBLEMS="${ENABLE_PROBLEMS}"
  -Dmdcraft_THREADS_TYPE="${THREADS_TYPE}"
  -Dmdcraft_ENABLE_BPNN="${ENABLE_BPNN}"
  -Dmdcraft_ENABLE_DeePMD="${ENABLE_DEEPMD}"
  -Dmdcraft_ENABLE_MLIP4="${ENABLE_MLIP4}"
  -Dmdcraft_FETCH_MISSING_DEPS="${FETCH_MISSING_DEPS}"
  -Dmdcraft_PYTHON_INSTALL_LAYOUT="${PYTHON_INSTALL_LAYOUT}"
)

if [[ "${ENABLE_MPI}" == "ON" ]]; then
  CONFIGURE_ARGS+=(
    -DMPI_C_COMPILER="${MPI_C_COMPILER}"
    -DMPI_CXX_COMPILER="${MPI_CXX_COMPILER}"
    -DMPIEXEC_EXECUTABLE="${MPIEXEC_BIN}"
  )
fi

"${RUN_ENV[@]}" "${CMAKE_BIN}" \
  "${GENERATOR_ARGS[@]}" \
  "${CONFIGURE_ARGS[@]}" \
  "${EXTRA_CMAKE_ARGS[@]}"

if [[ "${CONFIGURE_ONLY}" == "ON" ]]; then
  exit 0
fi

echo "Building and installing MDcraft"
BUILD_ARGS=(--build "${BUILD_DIR}" --parallel)
if [[ -n "${BUILD_PARALLELISM}" ]]; then
  BUILD_ARGS+=("${BUILD_PARALLELISM}")
fi
"${RUN_ENV[@]}" "${CMAKE_BIN}" "${BUILD_ARGS[@]}"
"${RUN_ENV[@]}" "${CMAKE_BIN}" --install "${BUILD_DIR}"

echo
echo "MDcraft is available to Python from:"
"${PYTHON_BIN}" - <<'PY'
import sysconfig
print(sysconfig.get_path("platlib"))
PY
