get_filename_component(LIB_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(LIB ${CMAKE_PROJECT_NAME}-${LIB_NAME})
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")

if(NOT MPI_FOUND)
find_package(MPI REQUIRED COMPONENTS C CXX)
endif()

add_library(${LIB} SHARED ${SOURCES})
target_link_libraries(${LIB} PUBLIC
    ${CMAKE_PROJECT_NAME}-data
    ${CMAKE_PROJECT_NAME}-lattice
    ${CMAKE_PROJECT_NAME}-tools
    MPI::MPI_C MPI::MPI_CXX
)

mdcraft_install_runtime_library(${LIB})
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${HEADERS}
  DESTINATION include/${CMAKE_PROJECT_NAME}/${LIB_NAME}
)
