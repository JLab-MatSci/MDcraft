get_filename_component(SUBLIB ${CMAKE_CURRENT_SOURCE_DIR} NAME)
add_library(_${PROJECT_NAME}_${LIB}_${SUBLIB} MODULE _${SUBLIB}.cxx)
target_link_libraries(
               _${PROJECT_NAME}_${LIB}_${SUBLIB} 
        PRIVATE pybind11::module
        PUBLIC  ${PROJECT_NAME}-${LIB}
)
set_target_properties(_${PROJECT_NAME}_${LIB}_${SUBLIB}
        PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                   SUFFIX "${PYTHON_MODULE_EXTENSION}"
)
install(TARGETS _${PROJECT_NAME}_${LIB}_${SUBLIB}
        DESTINATION ${mdcraft_PYTHON_PACKAGE_DIR}/${LIB}/${SUBLIB}
        COMPONENT python
        EXCLUDE_FROM_ALL
)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
        DESTINATION ${mdcraft_PYTHON_PACKAGE_DIR}/${LIB}/${SUBLIB}
        COMPONENT python
        EXCLUDE_FROM_ALL
)
