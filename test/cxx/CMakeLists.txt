set(MDCRAFT_CXX_TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/domain/test_domain.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/neibs/test_grid.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/neibs/test_verlet_list.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/solver/test_ljs_cut.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/solver/test_argon_ball_md.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/solver/test_aluminum_bulk_eam.cxx
)

add_executable(run_tests ${MDCRAFT_CXX_TEST_SOURCES})

find_package(Threads REQUIRED)

set_target_properties(run_tests PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
mdcraft_configure_runtime_executable(run_tests)

target_compile_definitions(run_tests PRIVATE
  MDCRAFT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
)

target_include_directories(run_tests PRIVATE
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
)

# target_compile_options(run_tests PRIVATE -fno-rtti)
target_link_options(run_tests PRIVATE -Wl,--no-as-needed)

target_link_libraries(run_tests PRIVATE
  GTest::gtest_main
  GTest::gtest
  Threads::Threads
  ${PROJECT_NAME}-data
  ${PROJECT_NAME}-io
  ${PROJECT_NAME}-lattice
  ${PROJECT_NAME}-neibs
  ${PROJECT_NAME}-solver
  ${PROJECT_NAME}-tools
)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
  target_link_libraries(run_tests PRIVATE stdc++fs)
endif()

gtest_discover_tests(run_tests)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/LMP_BALL_130k_Ar
  ${CMAKE_BINARY_DIR}/bin/LMP_BALL_130k_Ar
  COPYONLY
)

configure_file(
  ${PROJECT_SOURCE_DIR}/test/python/solver/Al2009.eam.alloy
  ${CMAKE_BINARY_DIR}/bin/Al2009.eam.alloy
  COPYONLY
)

install(TARGETS run_tests DESTINATION bin)

install(
  FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/LMP_BALL_130k_Ar
    ${PROJECT_SOURCE_DIR}/test/python/solver/Al2009.eam.alloy
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
)
