# Link scripts
list( APPEND test_files
    saber_compare.sh
    saber_cpplint.py
    saber_doc_overview.sh
    saber_fit_function.py
    saber_links.ksh
    saber_parallel.sh
    saber_plot.py
    saber_plot/bump_avg.py
    saber_plot/bump_contour_centered.py
    saber_plot/bump_contour_positive.py
    saber_plot/bump_diag.py
    saber_plot/bump_lct_cor.py
    saber_plot/bump_normality.py
    saber_plot/bump_randomization.py
    saber_plot/bump_sampling_grids.py
    saber_plot/bump_umf.py
    saber_plot/quench.py
    saber_set_ref.sh
    saber_setup.sh
    saber_tar_data.sh
    saber_tar_ref.sh
    saber_valgrind.sh
)

file( MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/saber_plot )
foreach( FILENAME ${test_files} )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                             ${PROJECT_SOURCE_DIR}/tools/${FILENAME}
                             ${CMAKE_BINARY_DIR}/bin/${FILENAME} )
endforeach()
