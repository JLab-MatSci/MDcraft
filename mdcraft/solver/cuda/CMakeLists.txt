get_filename_component(SUBLIB ${CMAKE_CURRENT_SOURCE_DIR} NAME)

if (mdcraft_ENABLE_CUDA_VARIANT)
  set(HEADERS
      single_var.h
  )
elseif(mdcraft_ENABLE_CUDA_CRTP)
  set(HEADERS
      pair_crtp.h
      single_crtp.h
      multi_crtp.h
  )
endif()

install(FILES ${HEADERS}
  DESTINATION include/${PROJECT_NAME}/${LIB}/${SUBLIB}
)