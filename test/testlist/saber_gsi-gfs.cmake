# dirac_gsi_gfs_global
saber_add_test( TARGET saber_dirac_gsi_gfs_global_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_gsi_gfs_global.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )
