get_filename_component(LIB ${CMAKE_CURRENT_SOURCE_DIR} NAME)
add_library(_${PROJECT_NAME}_${LIB} MODULE _${LIB}.cxx)
target_link_libraries(
         _${PROJECT_NAME}_${LIB}
  PRIVATE pybind11::module
   PUBLIC ${PROJECT_NAME}-${LIB}
)
mdcraft_configure_python_extension(_${PROJECT_NAME}_${LIB} "${LIB}")
mdcraft_install_python_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py "${LIB}")
