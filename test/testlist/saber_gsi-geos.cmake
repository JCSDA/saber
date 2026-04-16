# dirac_gsi_geos_global
saber_add_test( TARGET saber_dirac_gsi_geos_global_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_gsi_geos_global.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_gsi_geos_global_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_gsi_geos_global.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_gsi_geos_global_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_gsi_geos_global.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_gsi_geos_global_opt_1
saber_add_test( TARGET saber_dirac_gsi_geos_global_opt_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_gsi_geos_global_opt_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_gsi_geos_global_opt_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_gsi_geos_global_opt_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_gsi_geos_global_opt_1_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_gsi_geos_global_opt_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_gsi_geos_global_opt_2
saber_add_test( TARGET saber_dirac_gsi_geos_global_opt_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_gsi_geos_global_opt_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_gsi_geos_global_opt_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_gsi_geos_global_opt_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_gsi_geos_global_opt_2_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_gsi_geos_global_opt_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# dirac_gsi_geos_global_opt_3
saber_add_test( TARGET saber_dirac_gsi_geos_global_opt_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_gsi_geos_global_opt_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_gsi_geos_global_opt_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_gsi_geos_global_opt_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

saber_add_test( TARGET saber_dirac_gsi_geos_global_opt_3_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_gsi_geos_global_opt_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# optimization_bump_hdiag_gsi_1
saber_add_test( TARGET saber_optimization_bump_hdiag_gsi_1_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/optimization_bump_hdiag_gsi_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_gsi_geos_global_opt_1_1-1 )

saber_add_test( TARGET saber_optimization_bump_hdiag_gsi_1_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/optimization_bump_hdiag_gsi_1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_gsi_geos_global_opt_1_2-1 )

# optimization_bump_hdiag_gsi_2
saber_add_test( TARGET saber_optimization_bump_hdiag_gsi_2_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/optimization_bump_hdiag_gsi_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_gsi_geos_global_opt_2_1-1 )

saber_add_test( TARGET saber_optimization_bump_hdiag_gsi_2_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/optimization_bump_hdiag_gsi_2.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_gsi_geos_global_opt_2_2-1 )

# optimization_bump_hdiag_gsi_3
saber_add_test( TARGET saber_optimization_bump_hdiag_gsi_3_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/optimization_bump_hdiag_gsi_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_gsi_geos_global_opt_3_1-1 )

saber_add_test( TARGET saber_optimization_bump_hdiag_gsi_3_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/optimization_bump_hdiag_gsi_3.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_dirac_gsi_geos_global_opt_3_2-1 )
