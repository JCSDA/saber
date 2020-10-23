# Link scripts
list( APPEND test_files
    saber_compare.sh
    saber_cpplint.py
    saber_doc_overview.sh
    saber_links.ksh
    saber_parallel.sh
    saber_plot.py
    saber_plot/avg.py
    saber_plot/contour_centered.py
    saber_plot/contour_positive.py
    saber_plot/diag.py
    saber_plot/lct_cor.py
    saber_plot/normality.py
    saber_plot/randomization.py
    saber_plot/sampling_grids.py
    saber_plot/umf.py
    saber_set_ref.sh
    saber_setup.sh
    saber_tar_data.sh
    saber_tar_ref.sh
    saber_valgrind.sh
)

if( oops_qg_FOUND )
    list( APPEND test_files saber_qg_links.sh )
endif()

file( MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/saber_plot )
foreach( FILENAME ${test_files} )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                             ${PROJECT_SOURCE_DIR}/tools/${FILENAME}
                             ${CMAKE_BINARY_DIR}/bin/${FILENAME} )
endforeach()
