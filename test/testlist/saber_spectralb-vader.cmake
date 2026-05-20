message( STATUS "- SPECTRALB VADER" )

# RANDOMIZATION TESTS NEEDED FOR THIS CMAKE

# randomization_sqrtspectralb_2
saber_add_test( TARGET saber_randomization_sqrtspectralb_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_sqrtspectralb_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_sqrtspectralb_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# DIRAC TESTS

# dirac_spectralb_gauss_vader_1
saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_spectralb_gauss_vader_2
saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_spectralb_gauss_vader_3
saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_1_1-1 )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_1_2-1 )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_1_1-2 )

# dirac_spectralb_gauss_vader_4
saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_4_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_4_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_4_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_spectralb_gauss_vader_5
saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_5_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_1_1-1 )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_5_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_1_2-1 )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_5_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_1_1-2 )

# dirac_spectralb_gauss_vader_6
saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_6_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_6_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_6_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_spectralb_gauss_vader_7
saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_7_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_1_1-1 )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_7_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_1_2-1 )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_7_1-2
                MPI 1
                OMP 2
		LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_1_1-2 )


# NOTE: This test relies on Atlas 0.46.0, because earlier versions do not support
#       using ectrans on a split communicator.
if( ATLAS_SUPPORTS_SPECTRAL_ON_SPLIT_COMM )
  # dirac_spectral_parallel_gauss_vader_6
  # (should give similar results to dirac_spectral_gauss_vader_6_2_1)
  saber_add_test( TARGET saber_dirac_spectralb_parallel_gauss_vader_6_4-2
                  MPI 4
                  OMP 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_spectralb_parallel_gauss_vader_6.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x )

  # dirac_spectral_parallel_gauss_vader_7
  # (should give similar results to dirac_spectral_gauss_vader_7_1_1)
  saber_add_test( TARGET saber_dirac_spectralb_parallel_gauss_vader_7_7-1
                  MPI 7
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_subcomm_error_covariance_toolbox.x
                  ARGS testinput/dirac_spectralb_parallel_gauss_vader_7.yaml
                  DEPENDS saber_quench_subcomm_error_covariance_toolbox.x
                  TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_1_1-1 )
endif()

# dirac_write_fields
saber_add_test( TARGET saber_dirac_write_fields_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_write_fields.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_write_fields_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_write_fields.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_write_fields_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_write_fields.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# DIAGNOSTIC TESTS

# compare_diagnostics_gauss_vader
saber_add_test( TARGET saber_compare_diagnostics_gauss_vader
                TYPE SCRIPT
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_compare_dirac_diagnostics.py
                ARGS testinput/compare_diagnostics_gauss_vader.yaml )

# compare_diagnostics_outer_vars
saber_add_test( TARGET saber_compare_diagnostics_outer_vars
                TYPE SCRIPT
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_compare_dirac_diagnostics.py
                ARGS testinput/compare_diagnostics_outer_vars.yaml )

# ERROR COVARIANCE TRAINING MIO TESTS

# error_covariance_training_spectralb_mio
saber_add_test( TARGET saber_error_covariance_training_spectralb_mio_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_spectralb_mio.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_spectralb_mio_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_spectralb_mio.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_spectralb_mio_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_spectralb_mio.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_2_1-2 )
