# dirac_bifourier_ectrans
saber_add_test( TARGET saber_dirac_bifourier_ectrans_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_ectrans.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_ectrans_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_ectrans.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_ectrans_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_ectrans.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_bifourier_vordivtouv_ectrans_1
saber_add_test( TARGET saber_dirac_bifourier_vordivtouv_ectrans_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_vordivtouv_ectrans_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_vordivtouv_ectrans_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_vordivtouv_ectrans_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_vordivtouv_ectrans_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_vordivtouv_ectrans_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_bifourier_vordivtouv_ectrans_2
saber_add_test( TARGET saber_dirac_bifourier_vordivtouv_ectrans_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_vordivtouv_ectrans_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_vordivtouv_ectrans_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_vordivtouv_ectrans_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_vordivtouv_ectrans_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_vordivtouv_ectrans_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_bifourier_ectrans_1
saber_add_test( TARGET saber_error_covariance_training_bifourier_ectrans_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_ectrans_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_ectrans_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_ectrans_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_ectrans_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_ectrans_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_1-2 )
