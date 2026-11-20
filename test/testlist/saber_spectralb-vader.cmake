message( STATUS "- SPECTRALB VADER" )

# RANDOMIZATION TESTS NEEDED FOR THIS CMAKE

# randomization_sqrtspectralb_ens_analvars
saber_add_test( TARGET saber_randomization_sqrtspectralb_ens_analvars_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_ens_analvars.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_sqrtspectralb_ens_analvars_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_ens_analvars.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_sqrtspectralb_ens_analvars_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_ens_analvars.yaml
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

if (ATLAS_TEST_46_OR_GREATER)
  # dirac_spectralb_gauss_vader_with_smv_1
  saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_with_smv_1_1-1
                  MPI 1
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_spectralb_gauss_vader_with_smv_1.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x )

  saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_with_smv_1_2-1
                  MPI 2
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_spectralb_gauss_vader_with_smv_1.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x )

  saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_with_smv_1_1-2
                  MPI 1
                  OMP 2
                  LABELS   tier2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_spectralb_gauss_vader_with_smv_1.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x )
endif()

# dirac_spectralb_gauss_vader_3
saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

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
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_2_1-1 )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_7_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_2_2-1 )

saber_add_test( TARGET saber_dirac_spectralb_gauss_vader_7_1-2
                MPI 1
                OMP 2
		LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_vader_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_2_1-2 )

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

if(ATLAS_TEST_46_OR_GREATER)
  # compare_diagnostics_gauss_vader_with_smv
  saber_add_test( TARGET saber_compare_diagnostics_gauss_vader_with_smv
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_compare_dirac_diagnostics.py
                  ARGS testinput/compare_diagnostics_gauss_vader_with_smv.yaml )
endif()

# compare_diagnostics_gauss_vader
saber_add_test( TARGET saber_compare_diagnostics_gauss_vader
                TYPE SCRIPT
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_compare_dirac_diagnostics.py
                ARGS testinput/compare_diagnostics_gauss_vader.yaml )

# ERROR COVARIANCE TRAINING MIO TESTS

# error_covariance_training_spectralb_mio
saber_add_test( TARGET saber_error_covariance_training_spectralb_mio_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_spectralb_mio.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_ens_analvars_1-1 )

saber_add_test( TARGET saber_error_covariance_training_spectralb_mio_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_spectralb_mio.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_ens_analvars_2-1 )

saber_add_test( TARGET saber_error_covariance_training_spectralb_mio_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_spectralb_mio.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_ens_analvars_1-2 )

if(ATLAS_TEST_46_OR_GREATER)
  # randomization_sqrtspectralb_ens_analvars_with_smv
  saber_add_test( TARGET saber_randomization_sqrtspectralb_ens_analvars_with_smv_1-1
                  MPI 1
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/randomization_sqrtspectralb_ens_analvars_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x )

  saber_add_test( TARGET saber_randomization_sqrtspectralb_ens_analvars_with_smv_2-1
                  MPI 2
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/randomization_sqrtspectralb_ens_analvars_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x )

  saber_add_test( TARGET saber_randomization_sqrtspectralb_ens_analvars_with_smv_1-2
                  MPI 1
                  OMP 2
                  LABELS   tier2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/randomization_sqrtspectralb_ens_analvars_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x )

  # dirac_spectral_parallel_gauss_vader_with_smv_6
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

  # dirac_write_fields_with_smv
  saber_add_test( TARGET saber_dirac_write_fields_with_smv_1-1
                  MPI 1
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_write_fields_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x )

  saber_add_test( TARGET saber_dirac_write_fields_with_smv_2-1
                  MPI 2
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_write_fields_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x )

  saber_add_test( TARGET saber_dirac_write_fields_with_smv_1-2
                  MPI 1
                  OMP 2
                  LABELS   tier2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_write_fields_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x )

  # error_covariance_training_spectralb_mio_with_smv
  saber_add_test( TARGET saber_error_covariance_training_spectralb_mio_with_smv_1-1
                  MPI 1
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/error_covariance_training_spectralb_mio_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_randomization_sqrtspectralb_ens_analvars_with_smv_1-1 )

  saber_add_test( TARGET saber_error_covariance_training_spectralb_mio_with_smv_2-1
                  MPI 2
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/error_covariance_training_spectralb_mio_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_randomization_sqrtspectralb_ens_analvars_with_smv_2-1 )

  saber_add_test( TARGET saber_error_covariance_training_spectralb_mio_with_smv_1-2
                  MPI 1
                  OMP 2
                  LABELS   tier2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/error_covariance_training_spectralb_mio_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_randomization_sqrtspectralb_ens_analvars_with_smv_1-2 )
endif()
