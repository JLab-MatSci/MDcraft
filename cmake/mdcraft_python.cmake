function(mdcraft_require_full_install)
  install(CODE "
    if(NOT MDCRAFT_FULL_INSTALL)
      message(FATAL_ERROR \"Install MDcraft from the TOP-LEVEL build directory: cmake --install <build-root> (or make install there). Partial subdirectory installation can mix incompatible core and Python modules.\")
    endif()
  ")
endfunction()

function(mdcraft_require_python_module module_name)
  set(options)
  set(oneValueArgs PURPOSE INSTALL_HINT)
  cmake_parse_arguments(MDCRAFT_PYMOD "${options}" "${oneValueArgs}" "" ${ARGN})

  if(NOT MDCRAFT_PYMOD_PURPOSE)
    set(MDCRAFT_PYMOD_PURPOSE "this MDcraft configuration")
  endif()

  execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import importlib; importlib.import_module('${module_name}')"
    RESULT_VARIABLE _python_module_result
    OUTPUT_QUIET
    ERROR_VARIABLE _python_module_error
  )

  if(NOT _python_module_result EQUAL 0)
    string(REPLACE "\n" "\n  " _python_module_error "${_python_module_error}")
    set(_python_module_hint "")
    if(MDCRAFT_PYMOD_INSTALL_HINT)
      set(_python_module_hint "\nSuggested fix:\n  ${MDCRAFT_PYMOD_INSTALL_HINT}")
    endif()
    message(FATAL_ERROR
      "Python module '${module_name}' is required for ${MDCRAFT_PYMOD_PURPOSE}, "
      "but it is not available for '${Python_EXECUTABLE}'.${_python_module_hint}\n"
      "Import error:\n  ${_python_module_error}"
    )
  endif()
endfunction()

function(mdcraft_setup_python_environment)
  set(Python_FIND_VIRTUALENV FIRST)
  find_package(Python 3 REQUIRED COMPONENTS Interpreter Development.Module)
  set(Python_EXECUTABLE "${Python_EXECUTABLE}" CACHE FILEPATH "Python interpreter used to build mdcraft" FORCE)
  set(Python_EXECUTABLE "${Python_EXECUTABLE}" PARENT_SCOPE)

  mdcraft_require_python_module(
    numpy
    PURPOSE "the mdcraft Python bindings"
    INSTALL_HINT "\"${Python_EXECUTABLE}\" -m pip install numpy"
  )

  if(mdcraft_ENABLE_TESTS OR mdcraft_ENABLE_PROBLEMS)
    mdcraft_require_python_module(
     scipy
     PURPOSE "the installed Python tests and example scripts"
     INSTALL_HINT "\"${Python_EXECUTABLE}\" -m pip install scipy"
    )
    mdcraft_require_python_module(
     matplotlib
     PURPOSE "the installed Python tests and example scripts"
     INSTALL_HINT "\"${Python_EXECUTABLE}\" -m pip install matplotlib"
    )
    mdcraft_require_python_module(
     h5py
     PURPOSE "the installed Python tests and example scripts"
     INSTALL_HINT "\"${Python_EXECUTABLE}\" -m pip install h5py"
    )
    mdcraft_require_python_module(
     h5py
     PURPOSE "the installed Python tests and example scripts"
     INSTALL_HINT "\"${Python_EXECUTABLE}\" -m pip install ase"
    )
  endif()

  if(mdcraft_ENABLE_MPI)
    mdcraft_require_python_module(
      mpi4py
      PURPOSE "the MPI-enabled mdcraft Python bindings"
      INSTALL_HINT "\"${Python_EXECUTABLE}\" -m pip install mpi4py"
    )
    find_package(MPI4PY REQUIRED)
  endif()

  execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX') or '')"
    OUTPUT_VARIABLE _python_extension_suffix
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _python_extension_result
  )
  if(NOT _python_extension_result EQUAL 0 OR _python_extension_suffix STREQUAL "")
    set(_python_extension_suffix "${CMAKE_SHARED_MODULE_SUFFIX}")
  endif()

  execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_path('platlib'))"
    OUTPUT_VARIABLE _python_sitearch OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _python_sitearch_result)
  if(NOT _python_sitearch_result EQUAL 0 OR _python_sitearch STREQUAL "")
    message(FATAL_ERROR "Cannot determine site-packages for ${Python_EXECUTABLE}")
  endif()
  set(_python_package_root "${_python_sitearch}/${PROJECT_NAME}")
  set(_python_library_dir "${_python_package_root}/.libs")
  set(_python_extension_rpath "$ORIGIN/../.libs;$ORIGIN/../../.libs")

  set(PYTHON_MODULE_PREFIX "" CACHE INTERNAL "Prefix used for Python extension modules")
  set(PYTHON_MODULE_EXTENSION "${_python_extension_suffix}" CACHE INTERNAL "Suffix used for Python extension modules")
  set(mdcraft_PYTHON_PACKAGE_ROOT "${_python_package_root}" CACHE INTERNAL "Install root for the mdcraft Python package")
  set(mdcraft_PYTHON_LIBRARY_DIR "${_python_library_dir}" CACHE INTERNAL "Install directory for shared runtime libraries required by the Python package")
  set(mdcraft_PYTHON_EXTENSION_RPATH "${_python_extension_rpath}" CACHE INTERNAL "RPATH entries used by installed mdcraft Python extensions")
endfunction()

function(mdcraft_install_runtime_library target_name)
  mdcraft_require_full_install()
  set_target_properties(${target_name} PROPERTIES
    BUILD_RPATH_USE_ORIGIN ON
    INSTALL_RPATH "$ORIGIN"
    INSTALL_RPATH_USE_LINK_PATH TRUE
  )

  install(
    TARGETS ${target_name}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

  if(mdcraft_ENABLE_PYTHON)
    install(
      TARGETS ${target_name}
      LIBRARY DESTINATION "${mdcraft_PYTHON_LIBRARY_DIR}"
      ARCHIVE DESTINATION "${mdcraft_PYTHON_LIBRARY_DIR}"
      RUNTIME DESTINATION "${mdcraft_PYTHON_LIBRARY_DIR}"
    )
  endif()
endfunction()

function(mdcraft_configure_runtime_executable target_name)
  set_target_properties(${target_name} PROPERTIES
    BUILD_RPATH_USE_ORIGIN ON
    INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}"
    INSTALL_RPATH_USE_LINK_PATH TRUE
  )
endfunction()

function(mdcraft_configure_python_extension target_name relative_package_dir)
  mdcraft_require_full_install()
  string(REPLACE "/" "." _module_dir "${relative_package_dir}")
  set_property(GLOBAL APPEND PROPERTY MDCRAFT_INSTALL_MODULES
    "mdcraft.${_module_dir}.${target_name}")
  set_target_properties(${target_name} PROPERTIES
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
    BUILD_RPATH_USE_ORIGIN ON
    INSTALL_RPATH "${mdcraft_PYTHON_EXTENSION_RPATH}"
    INSTALL_RPATH_USE_LINK_PATH TRUE
  )

  install(
    TARGETS ${target_name}
    DESTINATION "${mdcraft_PYTHON_PACKAGE_ROOT}/${relative_package_dir}"
  )
endfunction()

function(mdcraft_install_python_file source_file relative_package_dir)
  mdcraft_require_full_install()
  install(
    FILES "${source_file}"
    DESTINATION "${mdcraft_PYTHON_PACKAGE_ROOT}/${relative_package_dir}"
  )
endfunction()

function(mdcraft_finalize_python_install)
  message(STATUS "MDcraft installation interpreter: ${Python_EXECUTABLE}")
  message(STATUS "MDcraft Python package destination: ${mdcraft_PYTHON_PACKAGE_ROOT}")
  message(STATUS "MDcraft Python core-library destination: ${mdcraft_PYTHON_LIBRARY_DIR}")
  message(STATUS "MDcraft C++ core-library destination: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
  get_property(_modules GLOBAL PROPERTY MDCRAFT_INSTALL_MODULES)
  set(MDCRAFT_INSTALL_MODULES "['mdcraft'")
  foreach(_module IN LISTS _modules)
    string(APPEND MDCRAFT_INSTALL_MODULES ", '${_module}'")
  endforeach()
  string(APPEND MDCRAFT_INSTALL_MODULES "]")
  configure_file("${PROJECT_SOURCE_DIR}/cmake/verify_python_install.py.in"
                 "${PROJECT_BINARY_DIR}/verify_python_install.py" @ONLY)
  if(NOT CMAKE_CROSSCOMPILING)
    install(CODE "
      execute_process(COMMAND \"${Python_EXECUTABLE}\" \"${PROJECT_BINARY_DIR}/verify_python_install.py\"
                      RESULT_VARIABLE _mdcraft_check)
      if(NOT _mdcraft_check EQUAL 0)
        message(FATAL_ERROR \"MDcraft Python installation verification failed. See the import/destination error above.\")
      endif()
    ")
  endif()
endfunction()
