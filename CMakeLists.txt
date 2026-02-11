cmake_minimum_required(VERSION 3.20)
project(mdcraft VERSION 1.0.0 LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Don't manually set flags, cmake will set them correctly
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Default build type at CMake configure time")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release)

# Selectable C++ standard (per-target via mdcraft::options)
set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard")
set_property(CACHE CMAKE_CXX_STANDARD PROPERTY STRINGS 11 14 17)

# Central "project options" target: if you update all the subdirs `CMakeLists.txt`,
# you can link to this target and get all the options at once
add_library(mdcraft_options INTERFACE)
add_library(mdcraft::options ALIAS mdcraft_options)

target_compile_definitions(mdcraft_options INTERFACE $<$<CONFIG:Debug>:DEBUG>)

# Options
option(mdcraft_ENABLE_MPI "Distributed computing with MPI" OFF)
option(mdcraft_ENABLE_PYTHON "Python support" OFF)
option(mdcraft_ENABLE_TESTS "Build tests" OFF)
option(mdcraft_ENABLE_PROBLEMS "Build problems" OFF)
option(mdcraft_ENABLE_BPNN "Behler-Parrinello NN potentials support" OFF)
option(mdcraft_ENABLE_DeePMD "DeePMD NN potentials support" OFF)
option(mdcraft_ENABLE_MLIP4 "MTP potentials support from MLIP4" OFF)

# chosse Multi-threading option
set(mdcraft_THREADS_TYPE "STD" CACHE STRING "Choose thread style: STD or TBB")
set_property(CACHE mdcraft_THREADS_TYPE PROPERTY STRINGS STD TBB)

if(mdcraft_THREADS_TYPE STREQUAL "TBB")
  message(STATUS "Use TBB")
  set(mdcraft_ENABLE_TBB ON)
  find_package(TBB REQUIRED)
  target_link_libraries(mdcraft_options INTERFACE TBB::tbb)
else()
  message(STATUS "Use standard threads")
  set(mdcraft_ENABLE_TBB OFF)
endif()

# Dependencies
find_package(Eigen3)

if(mdcraft_ENABLE_DeePMD)
  find_package(DeePMD REQUIRED)
endif()

if(mdcraft_ENABLE_MLIP4)
  find_package(mlip4 REQUIRED)
endif()

# Compatibility vars (if your subdirs still refer to these)
set(mdcraft_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(mdcraft_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(mdcraft_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(mdcraft_VERSION "${mdcraft_VERSION_MAJOR},${mdcraft_VERSION_MINOR},${mdcraft_VERSION_PATCH}")
add_subdirectory(mdcraft)

if(mdcraft_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

if(mdcraft_ENABLE_PROBLEMS)
  add_subdirectory(problems)
endif()

# this is where most of the change is
if(mdcraft_ENABLE_PYTHON)
  set(PYBIND11_FINDPYTHON ON)
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(pybind11 CONFIG REQUIRED)

  execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_path('platlib') or '')"
    OUTPUT_VARIABLE mdcraft_PYTHON_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  if(NOT mdcraft_PYTHON_SITE_PACKAGES)
    message(FATAL_ERROR "Failed to detect Python site-packages path from ${Python3_EXECUTABLE}")
  endif()

  set(mdcraft_PYTHON_INSTALL_DIR "${mdcraft_PYTHON_SITE_PACKAGES}" CACHE PATH
    "Python site-packages install directory")
  set(mdcraft_PYTHON_PACKAGE_DIR "${mdcraft_PYTHON_INSTALL_DIR}/${PROJECT_NAME}")

  add_subdirectory(python)
  add_custom_target(pythoninstall
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
    COMMAND "${CMAKE_COMMAND}" --install "${CMAKE_BINARY_DIR}" --component python --config $<CONFIG>
    COMMENT "Installing Python modules"
    USES_TERMINAL
  )
endif()
