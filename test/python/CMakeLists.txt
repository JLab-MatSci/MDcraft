if(mdcraft_ENABLE_PYTHON)
  add_test(
    NAME python-import-mdcraft
    COMMAND "${Python_EXECUTABLE}" -c "import mdcraft, mdcraft.data, mdcraft.lattice; print(mdcraft.__file__)"
  )

  add_test(
    NAME python-test-atoms
    COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/data/test-atoms.py"
  )
endif()

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DESTINATION test
        FILES_MATCHING
          PATTERN "*.py"
          PATTERN "*.alloy")

if(mdcraft_ENABLE_PYTHON)
  set(_eam_serial_arg)
  if(NOT mdcraft_ENABLE_MPI)
    set(_eam_serial_arg --serial)
  endif()
  add_test(NAME python-eam-analytic
    COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/problems/eam-regression/test-analytic.py"
            ${_eam_serial_arg} --threads 2 --steps 3)
  add_test(NAME python-eam-lifetime
    COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/problems/eam-regression/test-NiAl.py"
            ${_eam_serial_arg} --threads 2 --steps 3 --lifetime-regression
            --output "${PROJECT_BINARY_DIR}/eam-test-output")
endif()
