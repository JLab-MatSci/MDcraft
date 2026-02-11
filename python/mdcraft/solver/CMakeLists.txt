get_filename_component(LIB ${CMAKE_CURRENT_SOURCE_DIR} NAME)
add_library(_${PROJECT_NAME}_${LIB} MODULE _${LIB}.cxx)
target_link_libraries(
               _${PROJECT_NAME}_${LIB} 
        PRIVATE pybind11::module
        PUBLIC  ${PROJECT_NAME}-${LIB}
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

add_subdirectory(boundary)
add_subdirectory(potential)
add_subdirectory(stepper)
add_subdirectory(thermostat)
