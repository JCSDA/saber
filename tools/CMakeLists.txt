# Link scripts
list( APPEND test_files
    saber_cpplint.py
    saber_compare_dirac_diagnostics.py
    saber_doc_overview.sh
)

file( MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/saber_plot )
foreach( FILENAME ${test_files} )
    file( CREATE_LINK
          ${PROJECT_SOURCE_DIR}/tools/${FILENAME}
          ${CMAKE_BINARY_DIR}/bin/${FILENAME} SYMBOLIC )
endforeach()
