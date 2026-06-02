message( STATUS "- SPECTRALB" )

# Order here is set so that "ctest -V -R spectralb" will work

# ENSURING REFERENCES ARE CONSISTENT BETWEEN TESTS
# compare_diagnostics_spectralb_write_read
saber_add_test( TARGET saber_compare_diagnostics_spectralb_write_read
                TYPE SCRIPT
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_compare_dirac_diagnostics.py
                ARGS testinput/compare_diagnostics_spectralb_write_read.yaml )

# compare_diagnostics_spectralb_from_saber_file_1
saber_add_test( TARGET saber_compare_diagnostics_spectralb_from_saber_file_1
                TYPE SCRIPT
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_compare_dirac_diagnostics.py
                ARGS testinput/compare_diagnostics_spectralb_from_saber_file_1.yaml )

# compare_diagnostics_sqrtspectralb_1
saber_add_test( TARGET saber_compare_diagnostics_sqrtspectralb_1
                TYPE SCRIPT
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_compare_dirac_diagnostics.py
                ARGS testinput/compare_diagnostics_sqrtspectralb_1.yaml )

# compare_diagnostics_sqrtspectralb_2
saber_add_test( TARGET saber_compare_diagnostics_sqrtspectralb_2
                TYPE SCRIPT
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_compare_dirac_diagnostics.py
                ARGS testinput/compare_diagnostics_sqrtspectralb_2.yaml )

# ERROR COVARIANCE CONVERSION TEST               
# error_covariance_training_spectralb_1
saber_add_test( TARGET saber_error_covariance_training_spectralb_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_spectralb_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_spectralb_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_spectralb_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_spectralb_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_spectralb_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# RANDOMIZATION TESTS

# randomization_csdual_sqrtspectralb
saber_add_test( TARGET saber_randomization_csdual_sqrtspectralb_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_csdual_sqrtspectralb.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_csdual_sqrtspectralb_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_csdual_sqrtspectralb.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_csdual_sqrtspectralb_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_csdual_sqrtspectralb.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# randomization_sqrtspectralb_1
saber_add_test( TARGET saber_randomization_sqrtspectralb_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_sqrtspectralb_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_sqrtspectralb_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_sqrtspectralb_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_1-1 )

saber_add_test( TARGET saber_randomization_sqrtspectralb_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_2-1 )

saber_add_test( TARGET saber_randomization_sqrtspectralb_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_1-2 )

saber_add_test( TARGET saber_randomization_sqrtspectralb_4_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_1-1 )

saber_add_test( TARGET saber_randomization_sqrtspectralb_4_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_2-1 )

saber_add_test( TARGET saber_randomization_sqrtspectralb_4_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_1-2 )

# randomization_sqrtspectralb_5
saber_add_test( TARGET saber_randomization_sqrtspectralb_5_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_sqrtspectralb_5_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_sqrtspectralb_5_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_sqrtspectralb_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# STANDARD ERROR COVARIANCE TRAINING TESTS

# error_covariance_training_spectralb_2
saber_add_test( TARGET saber_error_covariance_training_spectralb_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_spectralb_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_3_1-1 )

saber_add_test( TARGET saber_error_covariance_training_spectralb_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_spectralb_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_3_2-1 )

saber_add_test( TARGET saber_error_covariance_training_spectralb_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_spectralb_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_3_1-2 )

# error_covariance_training_spectralb_3
saber_add_test( TARGET saber_error_covariance_training_spectralb_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_spectralb_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_4_1-1 )

saber_add_test( TARGET saber_error_covariance_training_spectralb_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_spectralb_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_4_2-1 )

saber_add_test( TARGET saber_error_covariance_training_spectralb_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_spectralb_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_4_1-2 )

# PROCESS PERTS TESTS

# process_perts_spectralb_from_csdual_states_1
saber_add_test( TARGET saber_process_perts_spectralb_from_csdual_states_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_csdual_states_1.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_csdual_sqrtspectralb_1-1 )

saber_add_test( TARGET saber_process_perts_spectralb_from_csdual_states_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_csdual_states_1.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_csdual_sqrtspectralb_2-1 )

saber_add_test( TARGET saber_process_perts_spectralb_from_csdual_states_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_csdual_states_1.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_csdual_sqrtspectralb_1-2 )

# process_perts_spectralb_from_csdual_states_2
saber_add_test( TARGET saber_process_perts_spectralb_from_csdual_states_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_csdual_states_2.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_csdual_sqrtspectralb_1-1 )

saber_add_test( TARGET saber_process_perts_spectralb_from_csdual_states_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_csdual_states_2.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_csdual_sqrtspectralb_2-1 )

saber_add_test( TARGET saber_process_perts_spectralb_from_csdual_states_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_csdual_states_2.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_csdual_sqrtspectralb_1-2 )

if(ATLAS_TEST_46_OR_GREATER)
  # process_perts_spectralb_from_csdual_states_1
  saber_add_test( TARGET saber_process_perts_spectralb_from_csdual_states_with_smv_1_1-1
                  MPI 1
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                  ARGS testinput/process_perts_spectralb_from_csdual_states_with_smv_1.yaml
                  DEPENDS saber_quench_process_perts.x
                  TEST_DEPENDS saber_randomization_csdual_sqrtspectralb_1-1 )

  saber_add_test( TARGET saber_process_perts_spectralb_from_csdual_states_with_smv_1_2-1
                  MPI 2
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                  ARGS testinput/process_perts_spectralb_from_csdual_states_with_smv_1.yaml
                  DEPENDS saber_quench_process_perts.x
                  TEST_DEPENDS saber_randomization_csdual_sqrtspectralb_2-1 )

  saber_add_test( TARGET saber_process_perts_spectralb_from_csdual_states_with_smv_1_1-2
                  MPI 1
                  OMP 2
                  LABELS   tier2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                  ARGS testinput/process_perts_spectralb_from_csdual_states_with_smv_1.yaml
                  DEPENDS saber_quench_process_perts.x
                  TEST_DEPENDS saber_randomization_csdual_sqrtspectralb_1-2 )

  # process_perts_spectralb_from_csdual_states_2
  saber_add_test( TARGET saber_process_perts_spectralb_from_csdual_states_with_smv_2_1-1
                  MPI 1
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                  ARGS testinput/process_perts_spectralb_from_csdual_states_with_smv_2.yaml
                  DEPENDS saber_quench_process_perts.x
                  TEST_DEPENDS saber_randomization_csdual_sqrtspectralb_1-1 )

  saber_add_test( TARGET saber_process_perts_spectralb_from_csdual_states_with_smv_2_2-1
                  MPI 2
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                  ARGS testinput/process_perts_spectralb_from_csdual_states_with_smv_2.yaml
                  DEPENDS saber_quench_process_perts.x
                  TEST_DEPENDS saber_randomization_csdual_sqrtspectralb_2-1 )

  saber_add_test( TARGET saber_process_perts_spectralb_from_csdual_states_with_smv_2_1-2
                  MPI 1
                  OMP 2
                  LABELS   tier2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                  ARGS testinput/process_perts_spectralb_from_csdual_states_with_smv_2.yaml
                  DEPENDS saber_quench_process_perts.x
                  TEST_DEPENDS saber_randomization_csdual_sqrtspectralb_1-2 )
endif()

# process_perts_spectralb_from_gauss_perts_1
saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_1.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_1_1-1 )

saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_1.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_1_2-1 )

saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_1.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_1_1-2 )

# process_perts_spectralb_from_gauss_perts_2
saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_2.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_1_1-1 )

saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_2.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_1_2-1 )

saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_2.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_1_1-2 )

# process_perts_spectralb_from_gauss_perts_3
saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_3.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_1_1-1 )

saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_3.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_1_2-1 )

saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_3.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_1_1-2 )

# process_perts_spectralb_from_gauss_perts_4
saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_4_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_4.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_5_1-1 )

saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_4_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_4.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_5_2-1 )

saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_4_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_4.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_5_1-2 )

# process_perts_spectralb_from_gauss_perts_5
saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_5_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_5.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_5_1-1 )

saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_5_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_5.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_5_2-1 )

saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_5_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_5.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_5_1-2 )

# process_perts_spectralb_from_gauss_perts_5
saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_6_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_6.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_4_1-1 )

saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_6_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_6.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_4_2-1 )

saber_add_test( TARGET saber_process_perts_spectralb_from_gauss_perts_6_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_spectralb_from_gauss_perts_6.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_sqrtspectralb_4_1-2 )

# DIRAC TESTS
if(ATLAS_TEST_46_OR_GREATER)
  # dirac_ens_both_geom_with_smv
  saber_add_test( TARGET saber_dirac_ens_both_geom_with_smv_1-1
                  MPI 1
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_ens_both_geom_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_with_smv_2_1-1 )

  saber_add_test( TARGET saber_dirac_ens_both_geom_with_smv_2-1
                  MPI 2
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_ens_both_geom_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_with_smv_2_2-1 )

  saber_add_test( TARGET saber_dirac_ens_both_geom_with_smv_1-2
                  MPI 1
                  OMP 2
                  LABELS   tier2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_ens_both_geom_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_with_smv_2_1-2 )

  # dirac_ens_model_geom
  saber_add_test( TARGET saber_dirac_ens_model_geom_with_smv_1-1
                  MPI 1
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_ens_model_geom_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_with_smv_2_1-1 )

  saber_add_test( TARGET saber_dirac_ens_model_geom_with_smv_2-1
                  MPI 2
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_ens_model_geom_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_with_smv_2_2-1 )

  saber_add_test( TARGET saber_dirac_ens_model_geom_with_smv_1-2
                  MPI 1
                  OMP 2
                  LABELS   tier2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_ens_model_geom_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_with_smv_2_1-2 )

  # dirac_ens_other_geom_with_smv_1
  saber_add_test( TARGET saber_dirac_ens_other_geom_with_smv_1_1-1
                  MPI 1
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_ens_other_geom_with_smv_1.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_with_smv_2_1-1 )

  saber_add_test( TARGET saber_dirac_ens_other_geom_with_smv_1_2-1
                  MPI 2
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_ens_other_geom_with_smv_1.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_with_smv_2_2-1 )

  saber_add_test( TARGET saber_dirac_ens_other_geom_with_smv_1_1-2
                  MPI 1
                  OMP 2
                  LABELS   tier2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_ens_other_geom_with_smv_1.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_with_smv_2_1-2 )

  # dirac_ens_other_geom_with_smv_2
  saber_add_test( TARGET saber_dirac_ens_other_geom_with_smv_2_1-1
                  MPI 1
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_ens_other_geom_with_smv_2.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_with_smv_2_1-1 )

  saber_add_test( TARGET saber_dirac_ens_other_geom_with_smv_2_2-1
                  MPI 2
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_ens_other_geom_with_smv_2.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_with_smv_2_2-1 )

  saber_add_test( TARGET saber_dirac_ens_other_geom_with_smv_2_1-2
                  MPI 1
                  OMP 2
                  LABELS   tier2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_ens_other_geom_with_smv_2.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_with_smv_2_1-2 )

  # This should give similar results to dirac_ens_other_geom_with_smv_2
  saber_add_test( TARGET saber_dirac_ens_parallel_other_geom_2_6-2
                  MPI 6
                  OMP 2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_ens_parallel_other_geom_2.yaml
		              LABELS   tier2
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_with_smv_2_1-2 )	

  # Compare model and other geom
  saber_add_test( TARGET saber_compare_diagnostics_ens_geom
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_compare_dirac_diagnostics.py
                  ARGS testinput/compare_diagnostics_ens_geom.yaml )

  saber_add_test( TARGET saber_compare_diagnostics_other_geom_parallel
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_compare_dirac_diagnostics.py
                  ARGS testinput/compare_diagnostics_other_geom_parallel.yaml)

  # dirac_interpolation_with_smv_1
  saber_add_test( TARGET saber_dirac_interpolation_with_smv_1_1-1
                  MPI 1
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_interpolation_with_smv_1.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x )

  saber_add_test( TARGET saber_dirac_interpolation_with_smv_1_2-1
                  MPI 2
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_interpolation_with_smv_1.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x )

  saber_add_test( TARGET saber_dirac_interpolation_with_smv_1_1-2
                  MPI 1
                  OMP 2
                  LABELS   tier2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_interpolation_with_smv_1.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x )


  # dirac_spectralb_from_CS_with_smv
  saber_add_test( TARGET saber_dirac_spectralb_from_CS_with_smv_1-1
                  MPI 1
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_spectralb_from_CS_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x )

  saber_add_test( TARGET saber_dirac_spectralb_from_CS_with_smv_2-1
                  MPI 2
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_spectralb_from_CS_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x )

  saber_add_test( TARGET saber_dirac_spectralb_from_CS_with_smv_1-2
                  MPI 1
                  OMP 2
                  LABELS   tier2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_spectralb_from_CS_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x )

  saber_add_test( TARGET saber_randomization_sqrtspectralb_with_smv_3_1-1
                  MPI 1
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/randomization_sqrtspectralb_with_smv_3.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_error_covariance_training_spectralb_1_1-1 )

  saber_add_test( TARGET saber_randomization_sqrtspectralb_with_smv_3_2-1
                  MPI 2
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/randomization_sqrtspectralb_with_smv_3.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_error_covariance_training_spectralb_1_2-1 )

  saber_add_test( TARGET saber_randomization_sqrtspectralb_with_smv_3_1-2
                  MPI 1
                  OMP 2
                  LABELS   tier2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/randomization_sqrtspectralb_with_smv_3.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_error_covariance_training_spectralb_1_1-2 )
endif()

# dirac_ens_both_geom
saber_add_test( TARGET saber_dirac_ens_both_geom_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_ens_both_geom.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_2_1-1 )

saber_add_test( TARGET saber_dirac_ens_both_geom_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_ens_both_geom.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_2_2-1 )

saber_add_test( TARGET saber_dirac_ens_both_geom_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_ens_both_geom.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_2_1-2 )

# dirac_ens_model_geom
saber_add_test( TARGET saber_dirac_ens_model_geom_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_ens_model_geom.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_2_1-1 )

saber_add_test( TARGET saber_dirac_ens_model_geom_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_ens_model_geom.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_2_2-1 )

saber_add_test( TARGET saber_dirac_ens_model_geom_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_ens_model_geom.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_2_1-2 )

# dirac_ens_other_geom_1
saber_add_test( TARGET saber_dirac_ens_other_geom_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_ens_other_geom_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_2_1-1 )

saber_add_test( TARGET saber_dirac_ens_other_geom_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_ens_other_geom_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_2_2-1 )

saber_add_test( TARGET saber_dirac_ens_other_geom_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_ens_other_geom_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_2_1-2 )

# dirac_ens_other_geom_2
saber_add_test( TARGET saber_dirac_ens_other_geom_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_ens_other_geom_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_2_1-1 )

saber_add_test( TARGET saber_dirac_ens_other_geom_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_ens_other_geom_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_2_2-1 )

saber_add_test( TARGET saber_dirac_ens_other_geom_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_ens_other_geom_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_csdual_states_2_1-2 )

# dirac_interpolation_1
saber_add_test( TARGET saber_dirac_interpolation_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_interpolation_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_interpolation_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_interpolation_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_interpolation_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_interpolation_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_parallel_hybrid_id
saber_add_test( TARGET saber_dirac_parallel_hybrid_id_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_parallel_hybrid_id.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_parallel_hybrid_id_4-1
                MPI 4
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_parallel_hybrid_id.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_parallel_hybrid_stddev
saber_add_test( TARGET saber_dirac_parallel_hybrid_stddev_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_parallel_hybrid_stddev.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_stddev_1_2-1
                             saber_error_covariance_training_stddev_2_2-1 )

# WARNING: dependency to 2-1 tests
saber_add_test( TARGET saber_dirac_parallel_hybrid_stddev_4-1
                MPI 4
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_parallel_hybrid_stddev.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_stddev_1_2-1
                             saber_error_covariance_training_stddev_2_2-1 )

# dirac_parallel_hybrid_ratio_stddev
saber_add_test( TARGET saber_dirac_parallel_hybrid_ratio_stddev_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_parallel_hybrid_ratio_stddev.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_stddev_1_2-1
                             saber_error_covariance_training_stddev_2_2-1 )

# WARNING: dependency to 2-1 tests
saber_add_test( TARGET saber_dirac_parallel_hybrid_ratio_stddev_4-1
                MPI 4
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_parallel_hybrid_ratio_stddev.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_stddev_1_2-1
                             saber_error_covariance_training_stddev_2_2-1 )

# WARNING: dependency to 2-1 tests
saber_add_test( TARGET saber_dirac_parallel_hybrid_ratio_stddev_5-1
                MPI 5
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_parallel_hybrid_ratio_stddev.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_stddev_1_2-1
                             saber_error_covariance_training_stddev_2_2-1 )

# dirac_spectralb
saber_add_test( TARGET saber_dirac_spectralb_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_spectralb_and_touv
saber_add_test( TARGET saber_dirac_spectralb_and_touv_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_and_touv.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_and_touv_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_and_touv.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_and_touv_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_and_touv.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_spectralb_correlation_profiles
saber_add_test( TARGET saber_dirac_spectralb_correlation_profiles_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_correlation_profiles.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_correlation_profiles_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_correlation_profiles.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_correlation_profiles_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_correlation_profiles.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_spectralb_correlation_rescaling
saber_add_test( TARGET saber_dirac_spectralb_correlation_rescaling_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_correlation_rescaling.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_spectralb_correlation_profiles_1-1 )

saber_add_test( TARGET saber_dirac_spectralb_correlation_rescaling_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_correlation_rescaling.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_spectralb_correlation_profiles_2-1 )

saber_add_test( TARGET saber_dirac_spectralb_correlation_rescaling_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_correlation_rescaling.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_spectralb_correlation_profiles_1-2 )

# dirac_spectralb_covariance_profiles
saber_add_test( TARGET saber_dirac_spectralb_covariance_profiles_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_covariance_profiles.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_covariance_profiles_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_covariance_profiles.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_covariance_profiles_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_covariance_profiles.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_spectralb_covariance_rescaling_1
saber_add_test( TARGET saber_dirac_spectralb_covariance_rescaling_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_covariance_rescaling_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_spectralb_covariance_profiles_1-1 )

saber_add_test( TARGET saber_dirac_spectralb_covariance_rescaling_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_covariance_rescaling_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_spectralb_covariance_profiles_2-1 )

saber_add_test( TARGET saber_dirac_spectralb_covariance_rescaling_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_covariance_rescaling_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_spectralb_covariance_profiles_1-2 )

# dirac_spectralb_covariance_rescaling_2
saber_add_test( TARGET saber_dirac_spectralb_covariance_rescaling_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_covariance_rescaling_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_spectralb_covariance_rescaling_1_1-1 )

saber_add_test( TARGET saber_dirac_spectralb_covariance_rescaling_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_covariance_rescaling_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_spectralb_covariance_rescaling_1_2-1 )

saber_add_test( TARGET saber_dirac_spectralb_covariance_rescaling_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_covariance_rescaling_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_spectralb_covariance_rescaling_1_1-2 )

if(ATLAS_TEST_46_OR_GREATER)
  # dirac_spectralb_correlation_rescaling_with_smv
  saber_add_test( TARGET saber_dirac_spectralb_correlation_rescaling_with_smv_1-1
                  MPI 1
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_spectralb_correlation_rescaling_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_dirac_spectralb_correlation_profiles_1-1 )

  saber_add_test( TARGET saber_dirac_spectralb_correlation_rescaling_with_smv_2-1
                  MPI 2
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_spectralb_correlation_rescaling_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_dirac_spectralb_correlation_profiles_2-1 )

  saber_add_test( TARGET saber_dirac_spectralb_correlation_rescaling_with_smv_1-2
                  MPI 1
                  OMP 2
                  LABELS   tier2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_spectralb_correlation_rescaling_with_smv.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_dirac_spectralb_correlation_profiles_1-2 )

  # dirac_spectralb_covariance_rescaling_1
  saber_add_test( TARGET saber_dirac_spectralb_covariance_rescaling_with_smv_1_1-1
                  MPI 1
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_spectralb_covariance_rescaling_with_smv_1.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_dirac_spectralb_covariance_profiles_1-1 )

  saber_add_test( TARGET saber_dirac_spectralb_covariance_rescaling_with_smv_1_2-1
                  MPI 2
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_spectralb_covariance_rescaling_with_smv_1.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_dirac_spectralb_covariance_profiles_2-1 )

  saber_add_test( TARGET saber_dirac_spectralb_covariance_rescaling_with_smv_1_1-2
                  MPI 1
                  OMP 2
                  LABELS   tier2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_spectralb_covariance_rescaling_with_smv_1.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_dirac_spectralb_covariance_profiles_1-2 )

  # dirac_spectralb_covariance_rescaling_2
  saber_add_test( TARGET saber_dirac_spectralb_covariance_rescaling_with_smv_2_1-1
                  MPI 1
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_spectralb_covariance_rescaling_with_smv_2.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_dirac_spectralb_covariance_rescaling_1_1-1 )

  saber_add_test( TARGET saber_dirac_spectralb_covariance_rescaling_with_smv_2_2-1
                  MPI 2
                  OMP 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_spectralb_covariance_rescaling_with_smv_2.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_dirac_spectralb_covariance_rescaling_1_2-1 )

  saber_add_test( TARGET saber_dirac_spectralb_covariance_rescaling_with_smv_2_1-2
                  MPI 1
                  OMP 2
                  LABELS   tier2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                  ARGS testinput/dirac_spectralb_covariance_rescaling_with_smv_2.yaml
                  DEPENDS saber_quench_error_covariance_toolbox.x
                  TEST_DEPENDS saber_dirac_spectralb_covariance_rescaling_1_1-2 )
endif()

# dirac_spectralb_duplicate_variables
saber_add_test( TARGET saber_dirac_spectralb_duplicate_variables_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_duplicate_variables.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_duplicate_variables_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_duplicate_variables.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_duplicate_variables_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_duplicate_variables.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_spectralb_from_CS
saber_add_test( TARGET saber_dirac_spectralb_from_CS_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_CS.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_from_CS_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_CS.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_from_CS_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_CS.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_spectralb_from_L15
saber_add_test( TARGET saber_dirac_spectralb_from_L15_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_L15.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_from_L15_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_L15.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_from_L15_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_L15.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_spectralb_from_saber_file_1
saber_add_test( TARGET saber_dirac_spectralb_from_saber_file_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_saber_file_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_1-1)

saber_add_test( TARGET saber_dirac_spectralb_from_saber_file_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_saber_file_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_2-1)              

saber_add_test( TARGET saber_dirac_spectralb_from_saber_file_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_saber_file_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x 
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_1-2)

saber_add_test( TARGET saber_dirac_spectralb_from_saber_file_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_saber_file_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_1-1 )

saber_add_test( TARGET saber_dirac_spectralb_from_saber_file_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_saber_file_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_2-1 )

saber_add_test( TARGET saber_dirac_spectralb_from_saber_file_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_saber_file_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_1-2 )

saber_add_test( TARGET saber_dirac_spectralb_from_saber_file_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_saber_file_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_1-1 )

saber_add_test( TARGET saber_dirac_spectralb_from_saber_file_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_saber_file_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_2-1 )

saber_add_test( TARGET saber_dirac_spectralb_from_saber_file_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_saber_file_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_1-2 )

saber_add_test( TARGET saber_dirac_spectralb_from_saber_file_4_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_saber_file_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_1-1 )

saber_add_test( TARGET saber_dirac_spectralb_from_saber_file_4_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_saber_file_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_2-1 )

saber_add_test( TARGET saber_dirac_spectralb_from_saber_file_4_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_from_saber_file_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_spectralb_1_1-2 )

# dirac_spectralb_gauss_1
saber_add_test( TARGET saber_dirac_spectralb_gauss_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_4_1-1 )

saber_add_test( TARGET saber_dirac_spectralb_gauss_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_4_2-1 )

saber_add_test( TARGET saber_dirac_spectralb_gauss_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_4_1-2 )

# dirac_spectralb_gauss_2
saber_add_test( TARGET saber_dirac_spectralb_gauss_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_5_1-1 )

saber_add_test( TARGET saber_dirac_spectralb_gauss_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_5_2-1 )

saber_add_test( TARGET saber_dirac_spectralb_gauss_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_gauss_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_process_perts_spectralb_from_gauss_perts_5_1-2 )

# dirac_spectralb_localization_1
saber_add_test( TARGET saber_dirac_spectralb_localization_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_localization_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_localization_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_localization_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_localization_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_localization_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_spectralb_localization_2
saber_add_test( TARGET saber_dirac_spectralb_localization_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_localization_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_localization_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_localization_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectralb_localization_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectralb_localization_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_spectraltouv
saber_add_test( TARGET saber_dirac_spectraltouv_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectraltouv.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectraltouv_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectraltouv.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_spectraltouv_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_spectraltouv.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_sqrtspectralb
saber_add_test( TARGET saber_dirac_sqrtspectralb_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_sqrtspectralb.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_sqrtspectralb_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_sqrtspectralb.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_sqrtspectralb_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_sqrtspectralb.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_sqrtspectralb_correl_1
saber_add_test( TARGET saber_dirac_sqrtspectralb_correl_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_sqrtspectralb_correl_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_sqrtspectralb_correl_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_sqrtspectralb_correl_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_sqrtspectralb_correl_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_sqrtspectralb_correl_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_vertproj
saber_add_test( TARGET saber_dirac_vertproj_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_vertproj.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_vertproj_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_vertproj.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_vertproj_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_vertproj.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )
