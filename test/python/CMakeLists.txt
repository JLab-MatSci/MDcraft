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
