list( APPEND interfaces_testinput
  testinput/interfaces_3densvar_bump.yaml
  testinput/interfaces_3dvar_bump.yaml
  testinput/interfaces_3dvar_hybrid_bump.yaml
  testinput/interfaces_4densvar_bump.yaml
  testinput/interfaces_4dvar_drplanczos_bump.yaml
  testinput/interfaces_4dvar_drplanczos_hybrid_bump.yaml
  testinput/interfaces_dirac_bump_cov.yaml
  testinput/interfaces_dirac_bump_hyb.yaml
  testinput/interfaces_dirac_bump_lct.yaml
  testinput/interfaces_dirac_bump_loc_3d.yaml
  testinput/interfaces_dirac_bump_loc_4d.yaml
  testinput/interfaces_parameters_bump_cov.yaml
  testinput/interfaces_parameters_bump_hyb.yaml
  testinput/interfaces_parameters_bump_lct.yaml
  testinput/interfaces_parameters_bump_loc_3d.yaml
  testinput/interfaces_parameters_bump_loc_4d.yaml
)

list( APPEND interfaces_testoutput
  testoutput/interfaces_3densvar_bump.test
  testoutput/interfaces_3dvar_bump.test
  testoutput/interfaces_3dvar_hybrid_bump.test
  testoutput/interfaces_4densvar_bump.test
  testoutput/interfaces_4dvar_drplanczos_bump.test
  testoutput/interfaces_4dvar_drplanczos_hybrid_bump.test
  testoutput/interfaces_dirac_bump_cov.test
  testoutput/interfaces_dirac_bump_hyb.test
  testoutput/interfaces_dirac_bump_lct.test
  testoutput/interfaces_dirac_bump_loc_3d.test
  testoutput/interfaces_dirac_bump_loc_4d.test
  testoutput/interfaces_parameters_bump_cov.test
  testoutput/interfaces_parameters_bump_hyb.test
  testoutput/interfaces_parameters_bump_lct.test
  testoutput/interfaces_parameters_bump_loc_3d.test
  testoutput/interfaces_parameters_bump_loc_4d.test
)

# Create data directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${interfaces_testinput})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create data directory for reference output and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
foreach(FILENAME ${interfaces_testoutput})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

ecbuild_add_resources( TARGET   interfaces_test_scripts
                       SOURCES_PACK
                       ${interfaces_testinput}
                     )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bump)

#####################################################################
# BUMP parameters tests
#####################################################################

ecbuild_add_test( TARGET test_interfaces_parameters_bump_cov_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/interfaces_estimate_parameters.x
                  ARGS testinput/parameters_bump_cov.yaml testoutput/parameters_bump_cov.log.out
                  DEPENDS interfaces_estimate_parameters.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run )

ecbuild_add_test( TARGET test_interfaces_parameters_bump_cov_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/parameters_bump_cov.log.out testoutput/parameters_bump_cov.test
                  TEST_DEPENDS test_interfaces_parameters_bump_cov_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_interfaces_parameters_bump_hyb_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/interfaces_estimate_parameters.x
                  ARGS testinput/parameters_bump_hyb.yaml testoutput/parameters_bump_hyb.log.out
                  DEPENDS interfaces_estimate_parameters.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run )

ecbuild_add_test( TARGET test_interfaces_parameters_bump_hyb_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/parameters_bump_hyb.log.out testoutput/parameters_bump_hyb.test
                  TEST_DEPENDS test_interfaces_parameters_bump_hyb_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_interfaces_parameters_bump_lct_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/interfaces_estimate_parameters.x
                  ARGS testinput/parameters_bump_lct.yaml testoutput/parameters_bump_lct.log.out
                  DEPENDS interfaces_estimate_parameters.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run )

ecbuild_add_test( TARGET test_interfaces_parameters_bump_lct_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/parameters_bump_lct.log.out testoutput/parameters_bump_lct.test
                  TEST_DEPENDS test_interfaces_parameters_bump_lct_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_interfaces_parameters_bump_loc_3d_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/interfaces_estimate_parameters.x
                  ARGS testinput/parameters_bump_loc_3d.yaml testoutput/parameters_bump_loc_3d.log.out
                  DEPENDS interfaces_estimate_parameters.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run )

ecbuild_add_test( TARGET test_interfaces_parameters_bump_loc_3d_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/parameters_bump_loc_3d.log.out testoutput/parameters_bump_loc_3d.test
                  TEST_DEPENDS test_interfaces_parameters_bump_loc_3d_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_interfaces_parameters_bump_loc_4d_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/interfaces_estimate_parameters.x
                  ARGS testinput/parameters_bump_loc_4d.yaml testoutput/parameters_bump_loc_4d.log.out
                  DEPENDS interfaces_estimate_parameters.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run )

ecbuild_add_test( TARGET test_interfaces_parameters_bump_loc_4d_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/parameters_bump_loc_4d.log.out testoutput/parameters_bump_loc_4d.test
                  TEST_DEPENDS test_interfaces_parameters_bump_loc_4d_run )

#####################################################################
# BUMP links
#####################################################################

ecbuild_add_test( TARGET test_interfaces_bump_links
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/bump_links.sh
                  TEST_DEPENDS test_interfaces_parameters_bump_cov_run
                               test_interfaces_parameters_bump_hyb_run
                               test_interfaces_parameters_bump_lct_run
                               test_interfaces_parameters_bump_loc_3d_run
                               test_interfaces_parameters_bump_loc_4d_run )

#####################################################################
# BUMP dirac tests
#####################################################################

ecbuild_add_test( TARGET test_interfaces_dirac_bump_cov_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_dirac.x
                  ARGS testinput/dirac_bump_cov.yaml testoutput/dirac_bump_cov.log.out
                  DEPENDS qg_dirac.x
                  TEST_DEPENDS test_qg_forecast_run test_interfaces_bump_links )

ecbuild_add_test( TARGET test_interfaces_dirac_bump_cov_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/dirac_bump_cov.log.out testoutput/dirac_bump_cov.test
                  TEST_DEPENDS test_interfaces_dirac_bump_cov_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_interfaces_dirac_bump_hyb_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_dirac.x
                  ARGS testinput/dirac_bump_hyb.yaml testoutput/dirac_bump_hyb.log.out
                  DEPENDS qg_dirac.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_interfaces_bump_links )

ecbuild_add_test( TARGET test_interfaces_dirac_bump_hyb_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/dirac_bump_hyb.log.out testoutput/dirac_bump_hyb.test
                  TEST_DEPENDS test_interfaces_dirac_bump_hyb_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_interfaces_dirac_bump_lct_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_dirac.x
                  ARGS testinput/dirac_bump_lct.yaml testoutput/dirac_bump_lct.log.out
                  DEPENDS qg_dirac.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_interfaces_bump_links )

ecbuild_add_test( TARGET test_interfaces_dirac_bump_lct_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/dirac_bump_lct.log.out testoutput/dirac_bump_lct.test
                  TEST_DEPENDS test_interfaces_dirac_bump_lct_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_interfaces_dirac_bump_loc_3d_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_dirac.x
                  ARGS testinput/dirac_bump_loc_3d.yaml testoutput/dirac_bump_loc_3d.log.out
                  DEPENDS qg_dirac.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_interfaces_bump_links )

ecbuild_add_test( TARGET test_interfaces_dirac_bump_loc_3d_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/dirac_bump_loc_3d.log.out testoutput/dirac_bump_loc_3d.test
                  TEST_DEPENDS test_interfaces_dirac_bump_loc_3d_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_interfaces_dirac_bump_loc_4d_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_dirac.x
                  ARGS testinput/dirac_bump_loc_4d.yaml testoutput/dirac_bump_loc_4d.log.out
                  DEPENDS qg_dirac.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_interfaces_bump_links )

ecbuild_add_test( TARGET test_interfaces_dirac_bump_loc_4d_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/dirac_bump_loc_4d.log.out testoutput/dirac_bump_loc_4d.test
                  TEST_DEPENDS test_interfaces_dirac_bump_loc_4d_run )

#####################################################################
# 3d variational tests
#####################################################################

ecbuild_add_test( TARGET test_interfaces_3densvar_bump_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/3densvar_bump.yaml testoutput/3densvar_bump.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_interfaces_make_obs_3d_run test_interfaces_bump_links )

ecbuild_add_test( TARGET test_interfaces_3densvar_bump_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/3densvar_bump.log.out testoutput/3densvar_bump.test
                  TEST_DEPENDS test_interfaces_3densvar_bump_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_interfaces_3dvar_bump_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/3dvar_bump.yaml testoutput/3dvar_bump.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_interfaces_make_obs_3d_run test_interfaces_bump_links )

ecbuild_add_test( TARGET test_interfaces_3dvar_bump_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/3dvar_bump.log.out testoutput/3dvar_bump.test
                  TEST_DEPENDS test_interfaces_3dvar_bump_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_interfaces_3dvar_hybrid_bump_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/3dvar_hybrid_bump.yaml testoutput/3dvar_hybrid_bump.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_interfaces_make_obs_3d_run test_interfaces_bump_links )

ecbuild_add_test( TARGET test_interfaces_3dvar_hybrid_bump_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/3dvar_hybrid_bump.log.out testoutput/3dvar_hybrid_bump.test
                  TEST_DEPENDS test_interfaces_3dvar_hybrid_bump_run )

#####################################################################
# 4d variational tests
#####################################################################

ecbuild_add_test( TARGET test_interfaces_4densvar_bump_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/4densvar_bump.yaml testoutput/4densvar_bump.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_interfaces_make_obs_4d_12h_run test_interfaces_bump_links )

ecbuild_add_test( TARGET test_interfaces_4densvar_bump_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4densvar_bump.log.out testoutput/4densvar_bump.test
                  TEST_DEPENDS test_interfaces_4densvar_bump_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_interfaces_4dvar_drplanczos_bump_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/4dvar_drplanczos_bump.yaml testoutput/4dvar_drplanczos_bump.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_interfaces_make_obs_4d_24h_run test_interfaces_bump_links )

ecbuild_add_test( TARGET test_interfaces_4dvar_drplanczos_bump_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar_drplanczos_bump.log.out testoutput/4dvar_drplanczos_bump.test
                  TEST_DEPENDS test_interfaces_4dvar_drplanczos_bump_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_interfaces_4dvar_drplanczos_hybrid_bump_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/4dvar_drplanczos_hybrid_bump.yaml testoutput/4dvar_drplanczos_hybrid_bump.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_interfaces_make_obs_4d_24h_run test_interfaces_bump_links )

ecbuild_add_test( TARGET test_interfaces_4dvar_drplanczos_hybrid_bump_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar_drplanczos_hybrid_bump.log.out testoutput/4dvar_drplanczos_hybrid_bump.test
                  TEST_DEPENDS test_interfaces_4dvar_drplanczos_hybrid_bump_run )
