# convertstate_constant
saber_add_test( TARGET saber_convertstate_constant
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_convertstate.x
                ARGS testinput/convertstate_constant.yaml
                DEPENDS saber_quench_convertstate.x )

# compare_diagnostics_ens_noloc
saber_add_test( TARGET saber_compare_diagnostics_ens_noloc
                TYPE SCRIPT
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_compare_dirac_diagnostics.py
                ARGS testinput/compare_diagnostics_ens_noloc.yaml )

# dirac_biperiodization_1
saber_add_test( TARGET saber_dirac_biperiodization_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_biperiodization_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_fastlam_1-1 )

saber_add_test( TARGET saber_dirac_biperiodization_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_biperiodization_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_fastlam_2-1 )

saber_add_test( TARGET saber_dirac_biperiodization_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_biperiodization_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_fastlam_1-2 )

# dirac_biperiodization_2
saber_add_test( TARGET saber_dirac_biperiodization_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_biperiodization_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_fastlam_1-1 )

saber_add_test( TARGET saber_dirac_biperiodization_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_biperiodization_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_fastlam_2-1 )

saber_add_test( TARGET saber_dirac_biperiodization_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_biperiodization_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_fastlam_1-2 )

# dirac_bump_1
saber_add_test( TARGET saber_dirac_bump_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertstate_constant
                             saber_error_covariance_training_bump_hdiag-nicas_1_1-1
                             saber_error_covariance_training_bump_stddev_3_1-1
                             saber_error_covariance_training_bump_vbal_1_1-1
                             saber_error_covariance_training_bump_wind_1_1-1
                             saber_randomization_bump_nicas_L10L2_1-1
                             saber_randomization_bump_nicas_L10L2T18_1-1 )

saber_add_test( TARGET saber_dirac_bump_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertstate_constant
                             saber_error_covariance_training_bump_hdiag-nicas_1_2-1
                             saber_error_covariance_training_bump_stddev_3_2-1
                             saber_error_covariance_training_bump_vbal_1_2-1
                             saber_error_covariance_training_bump_wind_1_2-1
                             saber_randomization_bump_nicas_L10L2_2-1
                             saber_randomization_bump_nicas_L10L2T18_2-1 )

saber_add_test( TARGET saber_dirac_bump_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertstate_constant
                             saber_error_covariance_training_bump_hdiag-nicas_1_1-2
                             saber_error_covariance_training_bump_stddev_3_1-2
                             saber_error_covariance_training_bump_vbal_1_1-2
                             saber_error_covariance_training_bump_wind_1_1-2
                             saber_randomization_bump_nicas_L10L2_1-2
                             saber_randomization_bump_nicas_L10L2T18_1-2 )

# dirac_bump_2
saber_add_test( TARGET saber_dirac_bump_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertstate_constant
                             saber_error_covariance_training_bump_hdiag-nicas_1_1-1
                             saber_error_covariance_training_bump_stddev_3_1-1
                             saber_error_covariance_training_bump_vbal_1_1-1
                             saber_error_covariance_training_bump_wind_1_1-1
                             saber_randomization_bump_nicas_L10L2_1-1
                             saber_randomization_bump_nicas_L10L2T18_1-1 )

saber_add_test( TARGET saber_dirac_bump_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertstate_constant
                             saber_error_covariance_training_bump_hdiag-nicas_1_2-1
                             saber_error_covariance_training_bump_stddev_3_2-1
                             saber_error_covariance_training_bump_vbal_1_2-1
                             saber_error_covariance_training_bump_wind_1_2-1
                             saber_randomization_bump_nicas_L10L2_2-1
                             saber_randomization_bump_nicas_L10L2T18_2-1 )

saber_add_test( TARGET saber_dirac_bump_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertstate_constant
                             saber_error_covariance_training_bump_hdiag-nicas_1_1-2
                             saber_error_covariance_training_bump_stddev_3_1-2
                             saber_error_covariance_training_bump_vbal_1_1-2
                             saber_error_covariance_training_bump_wind_1_1-2
                             saber_randomization_bump_nicas_L10L2_1-2
                             saber_randomization_bump_nicas_L10L2T18_1-2 )

# dirac_bump_3
saber_add_test( TARGET saber_dirac_bump_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertstate_constant
                             saber_error_covariance_training_bump_hdiag-nicas_1_1-1
                             saber_error_covariance_training_bump_stddev_3_1-1
                             saber_error_covariance_training_bump_vbal_1_1-1
                             saber_error_covariance_training_bump_wind_1_1-1
                             saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_dirac_bump_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertstate_constant
                             saber_error_covariance_training_bump_hdiag-nicas_1_2-1
                             saber_error_covariance_training_bump_stddev_3_2-1
                             saber_error_covariance_training_bump_vbal_1_2-1
                             saber_error_covariance_training_bump_wind_1_2-1
                             saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_dirac_bump_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertstate_constant
                             saber_error_covariance_training_bump_hdiag-nicas_1_1-2
                             saber_error_covariance_training_bump_stddev_3_1-2
                             saber_error_covariance_training_bump_vbal_1_1-2
                             saber_error_covariance_training_bump_wind_1_1-2
                             saber_randomization_bump_nicas_L10L2_1-2 )

# dirac_bump_4
saber_add_test( TARGET saber_dirac_bump_4_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertstate_constant
                             saber_error_covariance_training_bump_hdiag-nicas_1_1-1
                             saber_error_covariance_training_bump_stddev_3_1-1
                             saber_error_covariance_training_bump_vbal_1_1-1
                             saber_error_covariance_training_bump_wind_1_1-1
                             saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_dirac_bump_4_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertstate_constant
                             saber_error_covariance_training_bump_hdiag-nicas_1_2-1
                             saber_error_covariance_training_bump_stddev_3_2-1
                             saber_error_covariance_training_bump_vbal_1_2-1
                             saber_error_covariance_training_bump_wind_1_2-1
                             saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_dirac_bump_4_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_convertstate_constant
                             saber_error_covariance_training_bump_hdiag-nicas_1_1-2
                             saber_error_covariance_training_bump_stddev_3_1-2
                             saber_error_covariance_training_bump_vbal_1_1-2
                             saber_error_covariance_training_bump_wind_1_1-2
                             saber_randomization_bump_nicas_L10L2_1-2 )

# dirac_bump_5
saber_add_test( TARGET saber_dirac_bump_5_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1
                             saber_error_covariance_training_bump_hdiag-nicas_1_1-1
                             saber_error_covariance_training_bump_stddev_3_1-1 )

saber_add_test( TARGET saber_dirac_bump_5_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1
                             saber_error_covariance_training_bump_hdiag-nicas_1_2-1
                             saber_error_covariance_training_bump_stddev_3_2-1 )

saber_add_test( TARGET saber_dirac_bump_5_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2
                             saber_error_covariance_training_bump_hdiag-nicas_1_1-2
                             saber_error_covariance_training_bump_stddev_3_1-2 )

# dirac_bump_6
saber_add_test( TARGET saber_dirac_bump_6_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_dirac_bump_6_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_dirac_bump_6_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# dirac_bump_7
saber_add_test( TARGET saber_dirac_bump_7_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag-nicas_4_1-1
                             saber_error_covariance_training_bump_stddev_6_1-1
                             saber_error_covariance_training_bump_vbal_7_1-1 )

saber_add_test( TARGET saber_dirac_bump_7_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag-nicas_4_2-1
                             saber_error_covariance_training_bump_stddev_6_2-1
                             saber_error_covariance_training_bump_vbal_7_2-1 )

saber_add_test( TARGET saber_dirac_bump_7_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag-nicas_4_1-2
                             saber_error_covariance_training_bump_stddev_6_1-2
                             saber_error_covariance_training_bump_vbal_7_1-2 )

# dirac_bump_8
saber_add_test( TARGET saber_dirac_bump_8_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_8.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bump_8_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_8.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_bump_8_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_8.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_bump_9
saber_add_test( TARGET saber_dirac_bump_9_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_9.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag-nicas_1_1-1
                             saber_error_covariance_training_bump_wind_2_1-1 )

saber_add_test( TARGET saber_dirac_bump_9_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_9.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag-nicas_1_2-1
                             saber_error_covariance_training_bump_wind_2_2-1 )

saber_add_test( TARGET saber_dirac_bump_9_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_bump_9.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag-nicas_1_1-2
                             saber_error_covariance_training_bump_wind_2_1-2 )

# dirac_ens_noloc_4d
saber_add_test( TARGET saber_dirac_ens_noloc_4d_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_ens_noloc_4d.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1
                             saber_randomization_bump_nicas_L10L2T18_2-1 )

# dirac_localization_bump_1
saber_add_test( TARGET saber_dirac_localization_bump_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_localization_bump_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_dirac_localization_bump_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_localization_bump_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_dirac_localization_bump_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_localization_bump_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# dirac_localization_bump_2
saber_add_test( TARGET saber_dirac_localization_bump_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_localization_bump_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_dirac_localization_bump_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_localization_bump_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_dirac_localization_bump_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_localization_bump_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# dirac_localization_bump_3
saber_add_test( TARGET saber_dirac_localization_bump_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_localization_bump_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_dirac_localization_bump_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_localization_bump_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_dirac_localization_bump_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_localization_bump_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# dirac_localization_bump_4
saber_add_test( TARGET saber_dirac_localization_bump_4_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_localization_bump_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_dirac_localization_bump_4_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_localization_bump_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_dirac_localization_bump_4_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_localization_bump_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# dirac_localization_mixed_1
saber_add_test( TARGET saber_dirac_localization_mixed_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_localization_mixed_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_O10L2_1-1
                             saber_randomization_bump_nicas_F12L2_1-1 )

saber_add_test( TARGET saber_dirac_localization_mixed_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_localization_mixed_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_O10L2_2-1
                             saber_randomization_bump_nicas_F12L2_2-1 )

saber_add_test( TARGET saber_dirac_localization_mixed_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_localization_mixed_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_O10L2_1-2
                             saber_randomization_bump_nicas_F12L2_1-2 )

# dirac_oops_ens_noloc_4d
saber_add_test( TARGET saber_dirac_oops_ens_noloc_4d_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_oops_ens_noloc_4d.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1
                             saber_randomization_bump_nicas_L10L2T18_2-1 )

# dirac_shadowlevels_1
saber_add_test( TARGET saber_dirac_shadowlevels_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_shadowlevels_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_shadowlevels_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_shadowlevels_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_shadowlevels_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_shadowlevels_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_shadowlevels_2
saber_add_test( TARGET saber_dirac_shadowlevels_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_shadowlevels_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_shadowlevels_1_1-1 )

saber_add_test( TARGET saber_dirac_shadowlevels_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_shadowlevels_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_shadowlevels_1_2-1 )

saber_add_test( TARGET saber_dirac_shadowlevels_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_shadowlevels_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_shadowlevels_1_1-2 )

# dirac_shadowlevels_3
saber_add_test( TARGET saber_dirac_shadowlevels_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_shadowlevels_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_shadowlevels_1_1-1 )

saber_add_test( TARGET saber_dirac_shadowlevels_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_shadowlevels_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_shadowlevels_1_2-1 )

saber_add_test( TARGET saber_dirac_shadowlevels_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_shadowlevels_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_shadowlevels_1_1-2 )

# dirac_shadowlevels_4
saber_add_test( TARGET saber_dirac_shadowlevels_4_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_shadowlevels_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_shadowlevels_4_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_shadowlevels_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_shadowlevels_4_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_shadowlevels_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_shadowlevels_5
saber_add_test( TARGET saber_dirac_shadowlevels_5_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_shadowlevels_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_shadowlevels_4_1-1 )

saber_add_test( TARGET saber_dirac_shadowlevels_5_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_shadowlevels_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_shadowlevels_4_2-1 )

saber_add_test( TARGET saber_dirac_shadowlevels_5_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_shadowlevels_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_shadowlevels_4_1-2 )

# dirac_stddev_1
saber_add_test( TARGET saber_dirac_stddev_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_stddev_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_stddev_1_1-1 )

saber_add_test( TARGET saber_dirac_stddev_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_stddev_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_stddev_1_2-1 )

saber_add_test( TARGET saber_dirac_stddev_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_stddev_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_stddev_1_1-2 )

# dirac_stddev_2
saber_add_test( TARGET saber_dirac_stddev_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_stddev_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_stddev_2_1-1 )

saber_add_test( TARGET saber_dirac_stddev_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_stddev_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_stddev_2_2-1 )

saber_add_test( TARGET saber_dirac_stddev_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_stddev_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_stddev_2_1-2 )

# dirac_stddev_3
saber_add_test( TARGET saber_dirac_stddev_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_stddev_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_stddev_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_stddev_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_stddev_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_stddev_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_stddev_4
saber_add_test( TARGET saber_dirac_stddev_4_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_stddev_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_stddev_1_1-1 )

saber_add_test( TARGET saber_dirac_stddev_4_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_stddev_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_stddev_1_2-1 )

saber_add_test( TARGET saber_dirac_stddev_4_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_stddev_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_stddev_1_1-2 )

# error_covariance_training_bump_1
saber_add_test( TARGET saber_error_covariance_training_bump_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_bump_2
saber_add_test( TARGET saber_error_covariance_training_bump_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_bump_hdiag_1
saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_bump_hdiag_2
saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag_1_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag_1_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag_1_1-2 )

# error_covariance_training_bump_hdiag_3
saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_bump_hdiag_4
saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_4_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_4_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_4_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_bump_hdiag_5
saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_5_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L12L2_1-1
                             saber_randomization_bump_nicas_L10L2_1-1
                             saber_randomization_bump_nicas_L10L2_static_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_5_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L12L2_2-1
                             saber_randomization_bump_nicas_L10L2_2-1
                             saber_randomization_bump_nicas_L10L2_static_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_5_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L12L2_1-2
                             saber_randomization_bump_nicas_L10L2_1-2
                             saber_randomization_bump_nicas_L10L2_static_1-2 )

# error_covariance_training_bump_hdiag_6
saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_6_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L12L2_1-1
                             saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_6_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L12L2_2-1
                             saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_6_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L12L2_1-2
                             saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_bump_hdiag_7
saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_7_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L12L2_1-1
                             saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_7_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L12L2_2-1
                             saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_7_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L12L2_1-2
                             saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_bump_hdiag_8
saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_8_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_8.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L12L2_1-1
                             saber_randomization_bump_nicas_L10L2_static_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_8_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_8.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L12L2_2-1
                             saber_randomization_bump_nicas_L10L2_static_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_8_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_8.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L12L2_1-2
                             saber_randomization_bump_nicas_L10L2_static_1-2 )

# error_covariance_training_bump_hdiag_9
saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_9_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_9.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L12L2_1-1
                             saber_randomization_bump_nicas_L10L2_static_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_9_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_9.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L12L2_2-1
                             saber_randomization_bump_nicas_L10L2_static_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_9_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_9.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L12L2_1-2
                             saber_randomization_bump_nicas_L10L2_static_1-2 )

# error_covariance_training_bump_hdiag_10
saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_10_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_10.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L10_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_10_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_10.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L10_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_10_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_10.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L10_1-2 )

# error_covariance_training_bump_hdiag_11
saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_11_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_11.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_11_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_11.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_11_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_11.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_bump_hdiag_12
saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_12_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_12.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_12_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_12.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_12_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_12.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_bump_hdiag_13
saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_13_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_13.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_13_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_13.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag_13_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag_13.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_bump_hdiag-nicas_1
saber_add_test( TARGET saber_error_covariance_training_bump_hdiag-nicas_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag-nicas_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag-nicas_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag-nicas_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag-nicas_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag-nicas_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_bump_hdiag-nicas_2
saber_add_test( TARGET saber_error_covariance_training_bump_hdiag-nicas_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag-nicas_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag-nicas_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag-nicas_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag-nicas_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag-nicas_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_bump_hdiag-nicas_4
saber_add_test( TARGET saber_error_covariance_training_bump_hdiag-nicas_4_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag-nicas_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag-nicas_4_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag-nicas_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag-nicas_4_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag-nicas_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_bump_hdiag-nicas_5
saber_add_test( TARGET saber_error_covariance_training_bump_hdiag-nicas_5_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag-nicas_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_unstructured_global_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag-nicas_5_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag-nicas_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_unstructured_global_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_hdiag-nicas_5_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_hdiag-nicas_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_unstructured_global_1-2 )

# error_covariance_training_bump_nicas_1
saber_add_test( TARGET saber_error_covariance_training_bump_nicas_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_bump_nicas_2
saber_add_test( TARGET saber_error_covariance_training_bump_nicas_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag-nicas_1_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag-nicas_1_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag-nicas_1_1-2 )

# error_covariance_training_bump_nicas_3
saber_add_test( TARGET saber_error_covariance_training_bump_nicas_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag-nicas_1_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag-nicas_1_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag-nicas_1_1-2 )

# error_covariance_training_bump_nicas_4
saber_add_test( TARGET saber_error_covariance_training_bump_nicas_4_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag-nicas_1_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_4_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag-nicas_1_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_4_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag-nicas_1_1-2 )

# error_covariance_training_bump_nicas_5
saber_add_test( TARGET saber_error_covariance_training_bump_nicas_5_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_5_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_5_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_bump_nicas_6
saber_add_test( TARGET saber_error_covariance_training_bump_nicas_6_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag_1_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_6_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag_1_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_6_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_hdiag_1_1-2 )

# error_covariance_training_bump_nicas_7
saber_add_test( TARGET saber_error_covariance_training_bump_nicas_7_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_7_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_7_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_bump_nicas_8
saber_add_test( TARGET saber_error_covariance_training_bump_nicas_8_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_8.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_8_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_8.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_8_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_8.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_bump_nicas_9
saber_add_test( TARGET saber_error_covariance_training_bump_nicas_9_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_9.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_9_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_9.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_9_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_9.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_bump_nicas_10
saber_add_test( TARGET saber_error_covariance_training_bump_nicas_10_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_10.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_10_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_10.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_10_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_10.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_bump_nicas_11
saber_add_test( TARGET saber_error_covariance_training_bump_nicas_11_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_11.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_11_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_11.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_nicas_11_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_nicas_11.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_bump_stddev_1
saber_add_test( TARGET saber_error_covariance_training_bump_stddev_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_stddev_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_stddev_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_bump_stddev_2
saber_add_test( TARGET saber_error_covariance_training_bump_stddev_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_stddev_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_stddev_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_bump_stddev_3
saber_add_test( TARGET saber_error_covariance_training_bump_stddev_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_stddev_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_stddev_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_bump_stddev_4
saber_add_test( TARGET saber_error_covariance_training_bump_stddev_4_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_stddev_4_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_stddev_4_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_bump_stddev_5
saber_add_test( TARGET saber_error_covariance_training_bump_stddev_5_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_stddev_1_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_stddev_5_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_stddev_1_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_stddev_5_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_stddev_1_1-2 )

# error_covariance_training_bump_stddev_6
saber_add_test( TARGET saber_error_covariance_training_bump_stddev_6_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_stddev_6_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_stddev_6_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_stddev_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_bump_vbal_1
saber_add_test( TARGET saber_error_covariance_training_bump_vbal_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_vbal_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_vbal_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_bump_vbal_2
saber_add_test( TARGET saber_error_covariance_training_bump_vbal_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_vbal_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_vbal_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_bump_vbal_3
saber_add_test( TARGET saber_error_covariance_training_bump_vbal_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_vbal_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_vbal_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_bump_vbal_4
saber_add_test( TARGET saber_error_covariance_training_bump_vbal_4_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_vbal_3_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_vbal_4_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_vbal_3_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_vbal_4_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_4.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_vbal_3_1-2 )

# error_covariance_training_bump_vbal_5
saber_add_test( TARGET saber_error_covariance_training_bump_vbal_5_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_vbal_3_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_vbal_5_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_vbal_3_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_vbal_5_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_5.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_vbal_3_1-2 )

# error_covariance_training_bump_vbal_6
saber_add_test( TARGET saber_error_covariance_training_bump_vbal_6_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_vbal_3_1-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_vbal_6_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_vbal_3_2-1 )

saber_add_test( TARGET saber_error_covariance_training_bump_vbal_6_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_6.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_error_covariance_training_bump_vbal_3_1-2 )

# error_covariance_training_bump_vbal_7
saber_add_test( TARGET saber_error_covariance_training_bump_vbal_7_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_vbal_7_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_vbal_7_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_vbal_7.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_bump_wind_1
saber_add_test( TARGET saber_error_covariance_training_bump_wind_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_wind_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_wind_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_wind_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_wind_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_wind_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_bump_wind_2
saber_add_test( TARGET saber_error_covariance_training_bump_wind_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_wind_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_wind_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_wind_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_error_covariance_training_bump_wind_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_bump_wind_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# error_covariance_training_stddev_1
saber_add_test( TARGET saber_error_covariance_training_stddev_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_stddev_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_stddev_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_stddev_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_stddev_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_stddev_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# error_covariance_training_stddev_2
saber_add_test( TARGET saber_error_covariance_training_stddev_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_stddev_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-1 )

saber_add_test( TARGET saber_error_covariance_training_stddev_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_stddev_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_2-1 )

saber_add_test( TARGET saber_error_covariance_training_stddev_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/error_covariance_training_stddev_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_bump_nicas_L10L2_1-2 )

# process_perts_bump_nicas_1
saber_add_test( TARGET saber_process_perts_bump_nicas_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_bump_nicas_1.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_diffusion_2_1-1 )

saber_add_test( TARGET saber_process_perts_bump_nicas_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_bump_nicas_1.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_diffusion_2_2-1 )

saber_add_test( TARGET saber_process_perts_bump_nicas_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_bump_nicas_1.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_diffusion_2_1-2 )

# process_perts_bump_nicas_2
saber_add_test( TARGET saber_process_perts_bump_nicas_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_bump_nicas_2.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_diffusion_2_1-1
                             saber_process_perts_bump_nicas_1_1-1 )

saber_add_test( TARGET saber_process_perts_bump_nicas_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_bump_nicas_2.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_diffusion_2_2-1
                             saber_process_perts_bump_nicas_1_2-1 )

saber_add_test( TARGET saber_process_perts_bump_nicas_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_bump_nicas_2.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_diffusion_2_1-2
                             saber_process_perts_bump_nicas_1_1-2 )

# process_perts_bump_nicas_3
saber_add_test( TARGET saber_process_perts_bump_nicas_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_bump_nicas_3.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_diffusion_2_1-1 )

saber_add_test( TARGET saber_process_perts_bump_nicas_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_bump_nicas_3.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_diffusion_2_2-1 )

saber_add_test( TARGET saber_process_perts_bump_nicas_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_process_perts.x
                ARGS testinput/process_perts_bump_nicas_3.yaml
                DEPENDS saber_quench_process_perts.x
                TEST_DEPENDS saber_randomization_diffusion_2_1-2 )

# randomization_bump_nicas_F12L2
saber_add_test( TARGET saber_randomization_bump_nicas_F12L2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_F12L2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_F12L2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_F12L2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_F12L2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_F12L2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# randomization_bump_nicas_L10L2
saber_add_test( TARGET saber_randomization_bump_nicas_L10L2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_L10L2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_L10L2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_L10L2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_L10L2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_L10L2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# randomization_bump_nicas_L10L2T18
saber_add_test( TARGET saber_randomization_bump_nicas_L10L2T18_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_L10L2T18.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_L10L2T18_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_L10L2T18.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_L10L2T18_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_L10L2T18.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# randomization_bump_nicas_L10L2_static
saber_add_test( TARGET saber_randomization_bump_nicas_L10L2_static_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_L10L2_static.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_L10L2_static_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_L10L2_static.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_L10L2_static_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_L10L2_static.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# randomization_bump_nicas_L10L10
saber_add_test( TARGET saber_randomization_bump_nicas_L10L10_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_L10L10.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_L10L10_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_L10L10.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_L10L10_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_L10L10.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# randomization_bump_nicas_L12L2
saber_add_test( TARGET saber_randomization_bump_nicas_L12L2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_L12L2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_L12L2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_L12L2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_L12L2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_L12L2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# randomization_bump_nicas_unstructured_disc
saber_add_test( TARGET saber_randomization_bump_nicas_unstructured_disc_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_unstructured_disc.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_unstructured_disc_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_unstructured_disc.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_unstructured_disc_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_unstructured_disc.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# randomization_bump_nicas_unstructured_global
saber_add_test( TARGET saber_randomization_bump_nicas_unstructured_global_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_unstructured_global.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_unstructured_global_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_unstructured_global.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_unstructured_global_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_unstructured_global.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# randomization_bump_nicas_unstructured_rectangle
saber_add_test( TARGET saber_randomization_bump_nicas_unstructured_rectangle_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_unstructured_rectangle.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_unstructured_rectangle_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_unstructured_rectangle.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_randomization_bump_nicas_unstructured_rectangle_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_bump_nicas_unstructured_rectangle.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )
