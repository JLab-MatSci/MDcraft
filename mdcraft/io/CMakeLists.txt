get_filename_component(LIB ${CMAKE_CURRENT_SOURCE_DIR} NAME)
add_subdirectory(vtk)
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")

add_library(${PROJECT_NAME}-${LIB} SHARED
        ${vtk_SOURCES}
        ${SOURCES}
)

target_link_libraries(${PROJECT_NAME}-${LIB}
        PUBLIC ${PROJECT_NAME}-data
)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    target_link_libraries(${PROJECT_NAME}-${LIB} PUBLIC stdc++fs)
endif()

if (mdcraft_ENABLE_MPI)
    target_link_libraries(${PROJECT_NAME}-${LIB}
            PUBLIC ${PROJECT_NAME}-decomp)
endif()

mdcraft_install_runtime_library(${PROJECT_NAME}-${LIB})
