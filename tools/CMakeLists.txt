# Create Data directory for test input data
list( APPEND test_files
  cpplint.py
  saber_compare.sh
  saber_compare_specific.sh
  saber_links.ksh
  saber_oops_links.sh
  saber_set_ref.sh
  saber_tar_data.sh
  saber_tar_ref.sh
  saber_valgrind.sh
)

if (${OOPS_FOUND})
  list( APPEND test_files
    saber_oops_links.sh
  )
endif()

foreach(FILENAME ${test_files})
  execute_process( COMMAND ${CMAKE_COMMAND} -E copy
     ${PROJECT_SOURCE_DIR}/tools/${FILENAME}
     ${CMAKE_BINARY_DIR}/bin/${FILENAME} )
endforeach(FILENAME)
