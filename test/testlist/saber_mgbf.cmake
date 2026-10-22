message( STATUS "- MGBF" )

# dirac_mgbf_1
saber_add_test( TARGET saber_dirac_mgbf_1
                MPI 4
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_mgbf_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )
