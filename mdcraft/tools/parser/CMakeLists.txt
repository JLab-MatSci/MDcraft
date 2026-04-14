get_filename_component(LIB ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(SOURCES
      base.cxx
      maise_model.cxx
)

add_library(${PROJECT_NAME}-${LIB} SHARED ${SOURCES})

target_link_libraries(${PROJECT_NAME}-${LIB} 
  PUBLIC ${PROJECT_NAME}-data
         ${PROJECT_NAME}-tools
)

mdcraft_install_runtime_library(${PROJECT_NAME}-${LIB})

set(HEADERS
      base.h
      maise_model.h
)

install(FILES ${HEADERS}
  DESTINATION include/${PROJECT_NAME}/${LIB}
)
