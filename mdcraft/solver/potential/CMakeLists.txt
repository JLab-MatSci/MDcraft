get_filename_component(SUBLIB ${CMAKE_CURRENT_SOURCE_DIR} NAME)

if (mdcraft_ENABLE_CUDA)
  add_subdirectory(cuda)
endif()

file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
set(${SUBLIB}_SOURCES ${SOURCES} PARENT_SCOPE)
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${HEADERS}
  DESTINATION include/${PROJECT_NAME}/${LIB}/${SUBLIB}
)