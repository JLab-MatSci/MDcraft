find_package(MPI4PY REQUIRED)
get_filename_component(LIB ${CMAKE_CURRENT_SOURCE_DIR} NAME)
add_library(_${PROJECT_NAME}_${LIB} MODULE _${LIB}.cxx)
target_include_directories(
    _${PROJECT_NAME}_${LIB} 
    PRIVATE ${MPI4PY_INCLUDE_DIR}
)
if(NOT MPI_FOUND)
    find_package(MPI REQUIRED COMPONENTS C)
endif()
target_link_libraries(
         _${PROJECT_NAME}_${LIB}
  PRIVATE pybind11::module
   PUBLIC ${PROJECT_NAME}-${LIB}
    MPI::MPI_C
)
set_target_properties(_${PROJECT_NAME}_${LIB}
  PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
             SUFFIX "${PYTHON_MODULE_EXTENSION}"
)
install(TARGETS _${PROJECT_NAME}_${LIB}
        DESTINATION ${mdcraft_PYTHON_PACKAGE_DIR}/${LIB}
        COMPONENT python
        EXCLUDE_FROM_ALL
)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
        DESTINATION ${mdcraft_PYTHON_PACKAGE_DIR}/${LIB}
        COMPONENT python
        EXCLUDE_FROM_ALL
)
