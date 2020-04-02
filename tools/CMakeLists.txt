# Create Data directory for test input data
list( APPEND test_files
  bump_compare.sh
  bump_compare_specific.sh
  bump_links.ksh
  bump_parallel.sh
  bump_set_ref.sh
  bump_tar_data.sh
  bump_tar_ref.sh
  bump_valgrind.sh
  bump_wrapper.sh
)

foreach(FILENAME ${test_files})
  execute_process( COMMAND ${CMAKE_COMMAND} -E copy
     ${PROJECT_SOURCE_DIR}/tools/${FILENAME}
     ${CMAKE_BINARY_DIR}/bin/${FILENAME} )
endforeach(FILENAME)
