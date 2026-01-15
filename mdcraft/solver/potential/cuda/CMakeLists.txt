set(GPU_SOURCES
	LJs.cu
	EAM.cu
)

set(GPU_HEADERS
	LJs.cuh
	EAM.cuh
	base_potential.cuh
)
install(FILES ${GPU_HEADERS}
	DESTINATION include/${PROJECT_NAME}/${LIB}/${SUBLIB}/cuda
)

add_library(${PROJECT_NAME}-${LIB}-cuda SHARED
	${GPU_SOURCES}
)

target_link_libraries(${PROJECT_NAME}-${LIB}-cuda
	PRIVATE CUDA::cudart
)
install(TARGETS ${PROJECT_NAME}-${LIB}-cuda DESTINATION lib)