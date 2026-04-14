# устанавливаем в переменную имя библиотеки по имени папки
# write to a variable the library name by folder name
get_filename_component(LIB ${CMAKE_CURRENT_SOURCE_DIR} NAME)
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
add_library(${PROJECT_NAME}-${LIB} SHARED ${SOURCES})
target_include_directories(${PROJECT_NAME}-${LIB} PUBLIC
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
)
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${HEADERS}
  DESTINATION include/${PROJECT_NAME}/${LIB}
)
mdcraft_install_runtime_library(${PROJECT_NAME}-${LIB})

target_link_libraries(${PROJECT_NAME}-${LIB} PUBLIC ${PROJECT_NAME}-data)

# link for multi-threading
if (mdcraft_ENABLE_TBB)
    mdcraft_find_tbb_dependency()
    target_link_libraries(${PROJECT_NAME}-${LIB} PUBLIC TBB::tbb)
else ()
    find_package(Threads REQUIRED)
    target_link_libraries(${PROJECT_NAME}-${LIB} PUBLIC ${CMAKE_THREAD_LIBS_INIT})
endif ()

if (mdcraft_ENABLE_MPI)
    if(NOT MPI_FOUND)
        find_package(MPI REQUIRED COMPONENTS C CXX)
    endif()
    target_link_libraries(${PROJECT_NAME}-${LIB} PUBLIC MPI::MPI_C MPI::MPI_CXX)
endif()

if (mdcraft_ENABLE_BPNN)
    add_subdirectory(parser)
    add_subdirectory(bpnn)
endif()
