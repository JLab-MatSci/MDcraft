cmake_minimum_required(VERSION 3.18)
project(mdcraft VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_FLAGS_DEBUG      "-g")
set(CMAKE_CXX_FLAGS_RELEASE   "-O3")
# Select build type (Debug/Release) interactively
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Default build type at CMake configure time")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17 CACHE STRING "Default CXX standard at CMake configure time")
set_property(CACHE CMAKE_CXX_STANDARD PROPERTY STRINGS "11;14;17")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

set(mdcraft_VERSION_MAJOR         1)
set(mdcraft_VERSION_MINOR         0)
set(mdcraft_VERSION_PATCH         0)
set(mdcraft_VERSION "${mdcraft_VERSION_MAJOR},${mdcraft_VERSION_MINOR},${mdcraft_VERSION_PATCH}")

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/.local" CACHE PATH "Install prefix" FORCE)
endif()

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(GNUInstallDirs)

option(mdcraft_ENABLE_MPI            "Distributed computing with MPI"                ON)
option(mdcraft_ENABLE_PYTHON         "Python support"                                OFF)
option(mdcraft_ENABLE_TESTS          "Build tests"                                   OFF)
option(mdcraft_ENABLE_PROBLEMS       "Build problems"                                OFF)
option(mdcraft_ENABLE_BPNN           "Behler-Parrinello NN potentials support"       OFF)
option(mdcraft_ENABLE_DeePMD         "DeePmd NN potentials support"                  OFF)
option(mdcraft_ENABLE_MLIP4          "MTP potentials support from MLIP4"             OFF)
option(mdcraft_ENABLE_CUDA           "GPU CUDA acceleration"                         OFF)
option(mdcraft_FETCH_MISSING_DEPS    "Download missing third-party CMake dependencies automatically" ON)

set(mdcraft_PYTHON_INSTALL_LAYOUT "site-packages" CACHE STRING "Python install layout: site-packages or legacy")
set_property(CACHE mdcraft_PYTHON_INSTALL_LAYOUT PROPERTY STRINGS "site-packages;legacy")
set(mdcraft_PYTHON_SITEARCH "" CACHE PATH "Absolute Python site-packages directory to install the mdcraft package into")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions(-DDEBUG)
endif()

# chosse Multi-threading option
set(mdcraft_THREADS_TYPE "STD" CACHE STRING "Choose thread style: STD or TBB")
set_property(CACHE mdcraft_THREADS_TYPE PROPERTY STRINGS STD TBB)
if (${mdcraft_THREADS_TYPE} STREQUAL TBB)
  message(STATUS "Use TBB")
  set(mdcraft_ENABLE_TBB ON)
else ()
  message(STATUS "Use standard threads")
  set(mdcraft_ENABLE_TBB OFF)
endif ()

include(mdcraft_dependencies)
include(mdcraft_python)

mdcraft_find_eigen_dependency()

if(mdcraft_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS C CXX)
endif()

if (mdcraft_ENABLE_CUDA)
  option(mdcraft_ENABLE_CUDA_VARIANT
      "Use custom variant for CUDA device applications"     OFF)
  option(mdcraft_ENABLE_CUDA_CRTP
      "Use CRTP pattern for CUDA device applications"        ON)
  option(mdcraft_ENABLE_CUDA_THRUST
      "Use thrust library for CUDA device applications"      ON)

  if (mdcraft_ENABLE_CUDA_VARIANT AND mdcraft_ENABLE_CUDA_CRTP)
    message(FATAL_ERROR "Choose only one option for build with CUDA")
  endif()

  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)

  if (mdcraft_ENABLE_VARIANT_CUDA)
    include_directories(${PROJECT_NAME}/external)
  endif()

  # example for NVIDIA GeForce GTX 1080 (sm_61)
  # and for Tesla V100 (sm_70)
  set(CMAKE_CUDA_ARCHITECTURES 61 70)

  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)

  set(CMAKE_CUDA_FLAGS_DEBUG "-g -G")
  set(CMAKE_CUDA_FLAGS_RELEASE "-O3 --use_fast_math")

  message(STATUS "CUDA enabled. Architectures: ${CMAKE_CUDA_ARCHITECTURES}")

  if (mdcraft_ENABLE_CUDA_THRUST)
    find_package(Thrust REQUIRED)

    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda")
  endif()
endif()

if(mdcraft_ENABLE_DeePMD)
  find_package(DeePMD REQUIRED)
endif()

if(mdcraft_ENABLE_MLIP4)
  find_package(mlip4 REQUIRED)
endif()

add_subdirectory(mdcraft)

if(mdcraft_ENABLE_PYTHON)
  mdcraft_setup_python_environment()
  mdcraft_find_pybind11_dependency()
  add_subdirectory(python)
endif()

if(mdcraft_ENABLE_TESTS)
  mdcraft_find_gtest_dependency()
  include(GoogleTest)

  enable_testing()

  add_subdirectory(test)

endif()

if(mdcraft_ENABLE_PROBLEMS)
add_subdirectory(problems)
endif()
