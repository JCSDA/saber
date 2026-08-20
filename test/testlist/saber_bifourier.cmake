message( STATUS "- BIFOURIER" )

# convertcov_bifourier_balance_1
saber_add_test( TARGET saber_convertcov_bifourier_balance_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_balance_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_convertcov_bifourier_balance_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_balance_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_convertcov_bifourier_balance_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_balance_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# convertcov_bifourier_balance_2
saber_add_test( TARGET saber_convertcov_bifourier_balance_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_balance_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertcov_bifourier_balance_1_1-1 )

saber_add_test( TARGET saber_convertcov_bifourier_balance_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_balance_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertcov_bifourier_balance_1_2-1 )

saber_add_test( TARGET saber_convertcov_bifourier_balance_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_balance_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertcov_bifourier_balance_1_1-2 )

# convertcov_bifourier_balance_3
saber_add_test( TARGET saber_convertcov_bifourier_balance_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_balance_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertcov_bifourier_balance_2_1-1 )

saber_add_test( TARGET saber_convertcov_bifourier_balance_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_balance_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertcov_bifourier_balance_2_2-1 )

saber_add_test( TARGET saber_convertcov_bifourier_balance_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_balance_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertcov_bifourier_balance_2_1-2 )

# convertcov_bifourier_covariance_1
saber_add_test( TARGET saber_convertcov_bifourier_covariance_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_covariance_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_convertcov_bifourier_covariance_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_covariance_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_convertcov_bifourier_covariance_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_covariance_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# convertcov_bifourier_covariance_2
saber_add_test( TARGET saber_convertcov_bifourier_covariance_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_covariance_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertcov_bifourier_covariance_1_1-1 )

saber_add_test( TARGET saber_convertcov_bifourier_covariance_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_covariance_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertcov_bifourier_covariance_1_2-1 )

saber_add_test( TARGET saber_convertcov_bifourier_covariance_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_covariance_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertcov_bifourier_covariance_1_1-2 )

# convertcov_bifourier_covariance_3
saber_add_test( TARGET saber_convertcov_bifourier_covariance_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_covariance_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertcov_bifourier_covariance_2_1-1 )

saber_add_test( TARGET saber_convertcov_bifourier_covariance_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_covariance_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertcov_bifourier_covariance_2_2-1 )

saber_add_test( TARGET saber_convertcov_bifourier_covariance_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/convertcov_bifourier_covariance_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertcov_bifourier_covariance_2_1-2 )

# error_covariance_training_bifourier_covariance_1
saber_add_test( TARGET saber_error_covariance_training_bifourier_covariance_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_covariance_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_covariance_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_covariance_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_covariance_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_covariance_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_covariance_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_covariance_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_covariance_1-2 )

# error_covariance_training_bifourier_covariance_2
saber_add_test( TARGET saber_error_covariance_training_bifourier_covariance_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_covariance_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_covariance_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_covariance_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_covariance_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_covariance_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_covariance_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_covariance_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_covariance_1-2 )

# error_covariance_training_bifourier_1
saber_add_test( TARGET saber_error_covariance_training_bifourier_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_1-2 )

# error_covariance_training_bifourier_2
saber_add_test( TARGET saber_error_covariance_training_bifourier_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_1-2 )

# error_covariance_training_bifourier_3
saber_add_test( TARGET saber_error_covariance_training_bifourier_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_1-2 )

# error_covariance_training_bifourier_4
saber_add_test( TARGET saber_error_covariance_training_bifourier_4_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bifourier_3_1-1
                             saber_randomization_bifourier_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_4_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bifourier_3_2-1
                             saber_randomization_bifourier_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_4_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bifourier_3_1-2
                             saber_randomization_bifourier_1-2 )

# error_covariance_training_bifourier_5
saber_add_test( TARGET saber_error_covariance_training_bifourier_5_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_5_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_5_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_1-2 )

# error_covariance_training_bifourier_6
saber_add_test( TARGET saber_error_covariance_training_bifourier_6_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_6_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_6_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_1-2 )

# error_covariance_training_bifourier_7
saber_add_test( TARGET saber_error_covariance_training_bifourier_7_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_7_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bifourier_7_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_1-2 )

# error_covariance_training_bifourier_8
saber_add_test( TARGET saber_error_covariance_training_bifourier_8_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_8.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bifourier_8_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_8.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bifourier_8_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bifourier_8.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_bifourier_balance_1
saber_add_test( TARGET saber_dirac_bifourier_balance_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_balance_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertcov_bifourier_balance_3_1-1 )

saber_add_test( TARGET saber_dirac_bifourier_balance_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_balance_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertcov_bifourier_balance_3_2-1 )

saber_add_test( TARGET saber_dirac_bifourier_balance_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_balance_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertcov_bifourier_balance_3_1-2 )

# dirac_bifourier_balance_2
saber_add_test( TARGET saber_dirac_bifourier_balance_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_balance_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_balance_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_balance_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_balance_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_balance_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_bifourier_balance_3
saber_add_test( TARGET saber_dirac_bifourier_balance_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_balance_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_balance_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_balance_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_balance_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_balance_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_bifourier_covariance_1
saber_add_test( TARGET saber_dirac_bifourier_covariance_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_covariance_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_covariance_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_covariance_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_covariance_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_covariance_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_bifourier_covariance_2
saber_add_test( TARGET saber_dirac_bifourier_covariance_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_covariance_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_covariance_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_covariance_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_covariance_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_covariance_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_bifourier_covariance_3
saber_add_test( TARGET saber_dirac_bifourier_covariance_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_covariance_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_covariance_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_covariance_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_covariance_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_covariance_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_bifourier
saber_add_test( TARGET saber_dirac_bifourier_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_bifourier_ens
saber_add_test( TARGET saber_dirac_bifourier_ens_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_ens.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_1-1 )

saber_add_test( TARGET saber_dirac_bifourier_ens_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_ens.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_2-1 )

saber_add_test( TARGET saber_dirac_bifourier_ens_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_ens.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bifourier_1-2 )

# dirac_bifourier_gridtospectral
saber_add_test( TARGET saber_dirac_bifourier_gridtospectral_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_gridtospectral.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_gridtospectral_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_gridtospectral.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_gridtospectral_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_gridtospectral.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_bifourier_vordivtouv_1
saber_add_test( TARGET saber_dirac_bifourier_vordivtouv_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_vordivtouv_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_vordivtouv_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_vordivtouv_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_vordivtouv_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_vordivtouv_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_bifourier_vordivtouv_2
saber_add_test( TARGET saber_dirac_bifourier_vordivtouv_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_vordivtouv_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_vordivtouv_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_vordivtouv_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_vordivtouv_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_vordivtouv_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_bifourier_vordivtouv_3
saber_add_test( TARGET saber_dirac_bifourier_vordivtouv_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_vordivtouv_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_vordivtouv_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_vordivtouv_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bifourier_vordivtouv_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bifourier_vordivtouv_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# randomization_bifourier_covariance
saber_add_test( TARGET saber_randomization_bifourier_covariance_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bifourier_covariance.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bifourier_covariance_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bifourier_covariance.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bifourier_covariance_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bifourier_covariance.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# randomization_bifourier
saber_add_test( TARGET saber_randomization_bifourier_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bifourier.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bifourier_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bifourier.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bifourier_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bifourier.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )
