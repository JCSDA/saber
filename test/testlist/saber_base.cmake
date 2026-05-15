# dirac_diffusion_1
saber_add_test( TARGET saber_dirac_diffusion_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_diffusion_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_diffusion_1_1-1 )

saber_add_test( TARGET saber_dirac_diffusion_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_diffusion_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_diffusion_1_2-1 )

saber_add_test( TARGET saber_dirac_diffusion_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_diffusion_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_diffusion_1_1-2 )

# dirac_diffusion_2
saber_add_test( TARGET saber_dirac_diffusion_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_diffusion_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_diffusion_2_1-1 )

saber_add_test( TARGET saber_dirac_diffusion_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_diffusion_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_diffusion_2_2-1 )

saber_add_test( TARGET saber_dirac_diffusion_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_diffusion_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_diffusion_2_1-2 )

# dirac_diffusion_implicit_vt
# implicit-vertical diffusion is a per-column 1D Cholesky solve with no
# horizontal MPI communication, so a single MPI/OMP config exercises the
# new numerics; horizontal parallelism is covered by the explicit-vt tests.
saber_add_test( TARGET saber_dirac_diffusion_implicit_vt_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_diffusion_implicit_vt.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_diffusion_implicit_vt_1-1 )

# dirac_id
saber_add_test( TARGET saber_dirac_id_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_id.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_id_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_id.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_id_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_id.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_id_4d
saber_add_test( TARGET saber_dirac_id_4d_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_id_4d.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_id_4d_timecov
saber_add_test( TARGET saber_dirac_id_4d_timecov_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_id_4d_timecov.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_id_seq4d
saber_add_test( TARGET saber_dirac_id_seq4d_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_id_seq4d.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_id_seq4d_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_id_seq4d.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_id_seq4d_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_id_seq4d.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_id_seq4d_timecov
saber_add_test( TARGET saber_dirac_id_seq4d_timecov_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_id_seq4d_timecov.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_id_seq4d_timecov_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_id_seq4d_timecov.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_id_seq4d_timecov_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_id_seq4d_timecov.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_interpolation_2
saber_add_test( TARGET saber_dirac_interpolation_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_interpolation_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_interpolation_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_interpolation_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_interpolation_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_interpolation_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_interpolation_3
saber_add_test( TARGET saber_dirac_interpolation_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_interpolation_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_interpolation_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_interpolation_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_interpolation_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_interpolation_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_interpolation_4
saber_add_test( TARGET saber_dirac_interpolation_4_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_interpolation_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_interpolation_4_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_interpolation_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_interpolation_4_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_interpolation_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_interpolation_5
saber_add_test( TARGET saber_dirac_interpolation_5_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_interpolation_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_interpolation_5_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_interpolation_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_interpolation_5_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_interpolation_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_orographic_interp
saber_add_test( TARGET saber_dirac_orographic_interp_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_orographic_interp.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_orographic_interp_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_orographic_interp.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_orographic_interp_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_orographic_interp.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_duplicate_variables
saber_add_test( TARGET saber_dirac_duplicate_variables_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_duplicate_variables.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_duplicate_variables_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_duplicate_variables.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_duplicate_variables_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_duplicate_variables.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_diffusion_1
saber_add_test( TARGET saber_error_covariance_training_diffusion_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_diffusion_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_diffusion_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_diffusion_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_diffusion_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_diffusion_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_diffusion_2
saber_add_test( TARGET saber_error_covariance_training_diffusion_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_diffusion_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_diffusion_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_diffusion_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_diffusion_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_diffusion_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_diffusion_3
saber_add_test( TARGET saber_error_covariance_training_diffusion_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_diffusion_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_diffusion_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_diffusion_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_diffusion_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_diffusion_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_diffusion_implicit_vt
saber_add_test( TARGET saber_error_covariance_training_diffusion_implicit_vt_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_diffusion_implicit_vt.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# process_perts_diffusion_1
saber_add_test( TARGET saber_process_perts_diffusion_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_diffusion_1.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_diffusion_2_1-1
                             saber_error_covariance_training_diffusion_3_1-1 )

saber_add_test( TARGET saber_process_perts_diffusion_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_diffusion_1.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_diffusion_2_2-1
                             saber_error_covariance_training_diffusion_3_2-1 )

saber_add_test( TARGET saber_process_perts_diffusion_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_diffusion_1.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_diffusion_2_1-2
                             saber_error_covariance_training_diffusion_3_1-2 )

# randomization_diffusion_1
saber_add_test( TARGET saber_randomization_diffusion_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_diffusion_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_diffusion_1_1-1 )

saber_add_test( TARGET saber_randomization_diffusion_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_diffusion_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_diffusion_1_2-1 )

saber_add_test( TARGET saber_randomization_diffusion_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_diffusion_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_diffusion_1_1-2 )

# randomization_diffusion_2
saber_add_test( TARGET saber_randomization_diffusion_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_diffusion_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_diffusion_2_1-1 )

saber_add_test( TARGET saber_randomization_diffusion_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_diffusion_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_diffusion_2_2-1 )

saber_add_test( TARGET saber_randomization_diffusion_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_diffusion_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_diffusion_2_1-2 )

# randomization_increment_variables
saber_add_test( TARGET saber_randomization_increment_variables_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_increment_variables.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_increment_variables_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_increment_variables.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_increment_variables_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_increment_variables.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )
