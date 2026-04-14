get_filename_component(LIB ${CMAKE_CURRENT_SOURCE_DIR} NAME)
add_subdirectory(boundary)
add_subdirectory(potential)
add_subdirectory(stepper)
add_subdirectory(thermostat)
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")

add_library(${PROJECT_NAME}-${LIB} SHARED
  ${boundary_SOURCES}
  ${potential_SOURCES}
  ${stepper_SOURCES}
  ${thermostat_SOURCES}
  ${SOURCES}
)

target_link_libraries(${PROJECT_NAME}-${LIB}
  PUBLIC ${PROJECT_NAME}-data
         ${PROJECT_NAME}-lattice
         ${PROJECT_NAME}-neibs
         ${PROJECT_NAME}-tools
)

if (mdcraft_ENABLE_MPI)
    target_link_libraries(${PROJECT_NAME}-${LIB}
            PUBLIC ${PROJECT_NAME}-decomp)
endif()

if (mdcraft_ENABLE_BPNN)
    target_link_libraries(${PROJECT_NAME}-${LIB}
          PUBLIC ${PROJECT_NAME}-bpnn)
endif()

if (mdcraft_ENABLE_DeePMD)
target_link_libraries(${PROJECT_NAME}-${LIB}
  PUBLIC DeePMD::deepmd_c
)
endif()

if (mdcraft_ENABLE_MLIP4)
target_link_libraries(${PROJECT_NAME}-${LIB}
  PUBLIC mlip4::interface-lib
)
endif()

mdcraft_install_runtime_library(${PROJECT_NAME}-${LIB})

if (mdcraft_ENABLE_CUDA)
    add_subdirectory(cuda)
endif()
