# Variance accessor tests for saber block chains.

saber_add_test( TARGET saber_test_variance_id_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_test_variance.x
                ARGS testinput/test_variance_id.yaml
                DEPENDS saber_quench_test_variance.x )

saber_add_test( TARGET saber_test_variance_stddev_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_test_variance.x
                ARGS testinput/test_variance_stddev.yaml
                DEPENDS saber_quench_test_variance.x )

saber_add_test( TARGET saber_test_variance_diffusion_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_test_variance.x
                ARGS testinput/test_variance_diffusion.yaml
                DEPENDS saber_quench_test_variance.x
                TEST_DEPENDS error_covariance_training_diffusion_2_1-1 )

saber_add_test( TARGET saber_test_variance_diffusion_stddev_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_test_variance.x
                ARGS testinput/test_variance_diffusion_stddev.yaml
                DEPENDS saber_quench_test_variance.x
                TEST_DEPENDS error_covariance_training_diffusion_2_1-1 )

saber_add_test( TARGET saber_test_variance_ensemble_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_test_variance.x
                ARGS testinput/test_variance_ensemble.yaml
                DEPENDS saber_quench_test_variance.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_test_variance_bump_nicas_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_test_variance.x
                ARGS testinput/test_variance_bump_nicas.yaml
                DEPENDS saber_quench_test_variance.x )
