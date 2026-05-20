message( STATUS "- FASTLAM REGIONAL INTERPOLATION" )

# dirac_fastlam_10
saber_add_test( TARGET saber_dirac_fastlam_10_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_fastlam_10.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_fastlam_10_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_fastlam_10.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_fastlam_10_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_fastlam_10.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )
