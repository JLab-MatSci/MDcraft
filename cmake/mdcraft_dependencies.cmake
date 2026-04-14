include(FetchContent)

set(MDCRAFT_EIGEN_TAG "3.4.0" CACHE STRING "Pinned Eigen release used when mdcraft fetches Eigen automatically")
set(MDCRAFT_PYBIND11_TAG "v2.10.4" CACHE STRING "Pinned pybind11 release used when mdcraft fetches pybind11 automatically")
set(MDCRAFT_GTEST_TAG "release-1.12.1" CACHE STRING "Pinned googletest release used when mdcraft fetches it automatically")
set(MDCRAFT_TBB_TAG "v2021.13.0" CACHE STRING "Pinned oneTBB release used when mdcraft fetches TBB automatically")

function(mdcraft_require_fetch_mode dependency_name)
  if(mdcraft_FETCH_MISSING_DEPS)
    message(STATUS "${dependency_name} was not found locally. Fetching a pinned upstream release.")
  else()
    message(FATAL_ERROR "${dependency_name} was not found. Install it locally or re-run CMake with -Dmdcraft_FETCH_MISSING_DEPS=ON.")
  endif()
endfunction()

function(mdcraft_find_eigen_dependency)
  find_package(Eigen3 QUIET CONFIG)
  if(TARGET Eigen3::Eigen)
    return()
  endif()

  find_package(Eigen3 QUIET)
  if(TARGET Eigen3::Eigen)
    return()
  endif()

  mdcraft_require_fetch_mode("Eigen3")

  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
  set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)

  FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG ${MDCRAFT_EIGEN_TAG}
    GIT_SHALLOW TRUE
  )
  FetchContent_MakeAvailable(eigen)

  if(TARGET eigen AND NOT TARGET Eigen3::Eigen)
    add_library(Eigen3::Eigen ALIAS eigen)
  endif()
endfunction()

function(mdcraft_find_pybind11_dependency)
  set(PYBIND11_FINDPYTHON ON CACHE BOOL "" FORCE)

  execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import pybind11; print(pybind11.get_include())"
    OUTPUT_VARIABLE _pybind11_include_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _pybind11_include_result
  )
  if(_pybind11_include_result EQUAL 0 AND EXISTS "${_pybind11_include_dir}")
    add_library(pybind11::module INTERFACE IMPORTED)
    set_target_properties(pybind11::module PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${_pybind11_include_dir}"
      INTERFACE_LINK_LIBRARIES "Python::Module"
    )
    return()
  endif()

  find_package(pybind11 QUIET CONFIG)
  if(TARGET pybind11::module)
    return()
  endif()

  execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_VARIABLE _pybind11_cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _pybind11_result
  )
  if(_pybind11_result EQUAL 0 AND EXISTS "${_pybind11_cmake_dir}")
    set(pybind11_DIR "${_pybind11_cmake_dir}" CACHE PATH "pybind11 CMake package directory" FORCE)
    find_package(pybind11 QUIET CONFIG PATHS "${_pybind11_cmake_dir}" NO_DEFAULT_PATH)
    if(TARGET pybind11::module)
      return()
    endif()
  endif()

  mdcraft_require_fetch_mode("pybind11")

  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG ${MDCRAFT_PYBIND11_TAG}
    GIT_SHALLOW TRUE
  )
  FetchContent_MakeAvailable(pybind11)
endfunction()

function(mdcraft_find_gtest_dependency)
  if(mdcraft_ENABLE_MPI AND mdcraft_FETCH_MISSING_DEPS)
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest.git
      GIT_TAG ${MDCRAFT_GTEST_TAG}
      GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(googletest)
    return()
  endif()

  find_package(GTest QUIET CONFIG)
  if(TARGET GTest::gtest)
    return()
  endif()

  find_package(GTest QUIET)
  if(TARGET GTest::gtest)
    return()
  endif()

  mdcraft_require_fetch_mode("GTest")

  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG ${MDCRAFT_GTEST_TAG}
    GIT_SHALLOW TRUE
  )
  FetchContent_MakeAvailable(googletest)
endfunction()

function(mdcraft_find_tbb_dependency)
  find_package(TBB QUIET CONFIG)
  if(TARGET TBB::tbb)
    return()
  endif()

  find_package(TBB QUIET)
  if(TARGET TBB::tbb)
    return()
  endif()

  mdcraft_require_fetch_mode("TBB")

  set(TBB_TEST OFF CACHE BOOL "" FORCE)
  set(TBB_STRICT OFF CACHE BOOL "" FORCE)

  FetchContent_Declare(
    oneTBB
    GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
    GIT_TAG ${MDCRAFT_TBB_TAG}
    GIT_SHALLOW TRUE
  )
  FetchContent_MakeAvailable(oneTBB)

  if(NOT TARGET TBB::tbb)
    message(FATAL_ERROR "TBB was fetched, but target TBB::tbb is still unavailable.")
  endif()
endfunction()
