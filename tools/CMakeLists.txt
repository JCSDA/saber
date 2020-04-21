# Create Data directory for test input data
list( APPEND test_files
    bump_compare.sh
    bump_links.ksh
    bump_parallel.sh
    bump_valgrind.sh
    bump_wrapper.sh
    cpplint.py
    saber_set_ref.sh
    saber_tar_data.sh
    saber_tar_ref.sh
)

if (${OOPS_FOUND} AND ${ENABLE_OOPS_TOYMODELS})
   list( APPEND test_files qg_links.sh )
endif()

foreach(FILENAME ${test_files})
    execute_process( COMMAND ${CMAKE_COMMAND} -E copy
                             ${PROJECT_SOURCE_DIR}/tools/${FILENAME}
                             ${CMAKE_BINARY_DIR}/bin/${FILENAME} )
endforeach(FILENAME)
