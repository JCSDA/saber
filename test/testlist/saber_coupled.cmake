# coupled_dirac_id
saber_add_test( TARGET coupled_dirac_id_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_coupled_error_covariance_toolbox.x
                ARGS testinput/coupled_dirac_id.yaml
                DEPENDS saber_quench_coupled_error_covariance_toolbox.x
                TEST_DEPENDS error_covariance_training_diffusion_2 )

saber_add_test( TARGET coupled_dirac_id_1-2
                MPI 1
                OMP 2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_coupled_error_covariance_toolbox.x
                ARGS testinput/coupled_dirac_id.yaml
                DEPENDS saber_quench_coupled_error_covariance_toolbox.x
                TEST_DEPENDS error_covariance_training_diffusion_2 )
