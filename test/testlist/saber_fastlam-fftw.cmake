# dirac_fastlam-fftw_1
saber_add_test( TARGET saber_dirac_fastlam-fftw_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_fastlam-fftw_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertstate_lam
                             saber_randomization_bump_nicas_lam_1_1-1
                             saber_randomization_bump_nicas_lam_2_1-1 )

saber_add_test( TARGET saber_dirac_fastlam-fftw_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_fastlam-fftw_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertstate_lam
                             saber_randomization_bump_nicas_lam_1_2-1
                             saber_randomization_bump_nicas_lam_2_2-1 )

saber_add_test( TARGET saber_dirac_fastlam-fftw_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_fastlam-fftw_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertstate_lam
                             saber_randomization_bump_nicas_lam_1_1-2
                             saber_randomization_bump_nicas_lam_2_1-2 )

# dirac_fastlam-fftw_2
saber_add_test( TARGET saber_dirac_fastlam-fftw_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_fastlam-fftw_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_fastlam-fftw_1_1-1 )

saber_add_test( TARGET saber_dirac_fastlam-fftw_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_fastlam-fftw_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_fastlam-fftw_1_2-1 )

saber_add_test( TARGET saber_dirac_fastlam-fftw_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_fastlam-fftw_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_fastlam-fftw_1_1-2 )
