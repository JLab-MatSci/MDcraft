get_filename_component(SUBLIB ${CMAKE_CURRENT_SOURCE_DIR} NAME)

if (mdcraft_ENABLE_CUDA)
set_source_files_properties(_potential.cxx PROPERTIES LANGUAGE CUDA)
endif()

add_library(_${PROJECT_NAME}_${LIB}_${SUBLIB} MODULE _${SUBLIB}.cxx)
target_link_libraries(
               _${PROJECT_NAME}_${LIB}_${SUBLIB} 
        PRIVATE pybind11::module
        PUBLIC  ${PROJECT_NAME}-${LIB}
)
mdcraft_configure_python_extension(_${PROJECT_NAME}_${LIB}_${SUBLIB} "${LIB}/${SUBLIB}")
mdcraft_install_python_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py "${LIB}/${SUBLIB}")
