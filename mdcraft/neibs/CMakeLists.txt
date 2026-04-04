# устанавливаем в переменную имя библиотеки по имени папки
# write to a variable the library name according to folder name
get_filename_component(LIB ${CMAKE_CURRENT_SOURCE_DIR} NAME)
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
add_library(${PROJECT_NAME}-${LIB} SHARED ${SOURCES})
target_link_libraries(${PROJECT_NAME}-${LIB} 
  PUBLIC ${PROJECT_NAME}-data
         ${PROJECT_NAME}-lattice
         ${PROJECT_NAME}-tools
)
mdcraft_install_runtime_library(${PROJECT_NAME}-${LIB})
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${HEADERS}
  DESTINATION include/${PROJECT_NAME}/${LIB}
)
