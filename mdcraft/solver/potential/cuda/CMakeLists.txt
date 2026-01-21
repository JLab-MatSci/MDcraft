set(GPU_HEADERS
	LJs.cuh
	EAM.cuh
	base_potential.cuh
)
install(FILES ${GPU_HEADERS}
	DESTINATION include/${PROJECT_NAME}/${LIB}/${SUBLIB}/cuda
)

add_library(${PROJECT_NAME}-${LIB}-cuda INTERFACE)

target_include_directories(${PROJECT_NAME}-${LIB}-cuda INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

target_link_libraries(${PROJECT_NAME}-${LIB}-cuda
	INTERFACE CUDA::cudart
)

install(TARGETS ${PROJECT_NAME}-${LIB}-cuda DESTINATION lib)