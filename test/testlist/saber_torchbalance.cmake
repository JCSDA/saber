message( STATUS "- TORCHBALANCE" )

# dirac_torchbalance
saber_add_test( TARGET saber_dirac_torchbalance_1-1
                MPI 1
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_torchbalance.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS generate_torchbalance_test_emulator )

saber_add_test( TARGET saber_dirac_torchbalance_2-1
                MPI 2
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_torchbalance.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS generate_torchbalance_test_emulator )

saber_add_test( TARGET saber_dirac_torchbalance_1-2
                MPI 1
                OMP 2
                LABELS   tier2
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_torchbalance.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS generate_torchbalance_test_emulator )
