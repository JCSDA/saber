# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list(APPEND saber_test
  bump_cortrack
  bump_hdiag_cor_diag_all
  bump_hdiag_draw_type
  bump_hdiag_grid
  bump_hdiag_hyb-avg_common
  bump_hdiag_hyb-rnd_common
  bump_hdiag_ldwv
  bump_hdiag_loc_gau_approx
  bump_hdiag_loc_histograms
  bump_hdiag_mask_lat
  bump_hdiag_mask_stddev_lower
  bump_hdiag_mask_stddev_ncontig
  bump_hdiag_mask_stddev_upper
  bump_hdiag-nicas_cor_common_univariate
  bump_hdiag-nicas_cor_specific_univariate
  bump_hdiag-nicas_double_fit
  bump_hdiag-nicas_lhom
  bump_hdiag-nicas_loc_adv
  bump_hdiag-nicas_loc_adv_inv
  bump_hdiag-nicas_local_diag_var_filter
  bump_hdiag-nicas_loc_common
  bump_hdiag-nicas_loc_common_univariate
  bump_hdiag-nicas_loc_common_weighted
  bump_hdiag-nicas_loc_specific_multivariate
  bump_hdiag-nicas_loc_specific_univariate
  bump_hdiag-nicas_mask_check
  bump_hdiag-nicas_network
  bump_hdiag-nicas_nonunit_diag
  bump_hdiag-nicas_nrep
  bump_hdiag-nicas_rvflt
  bump_hdiag_optimality
  bump_lct-nicas_diagonal
  bump_lct-nicas_mask_check
  bump_lct-nicas_network
  bump_lct-nicas_one_scale
  bump_lct_two_scales
  bump_nicas_consistency
  bump_nicas_fast_sampling
  bump_nicas_mpicom_lsqrt_a
  bump_nicas_mpicom_lsqrt_b
  bump_nicas_mpicom_lsqrt_c
  bump_nicas_pos_def
  bump_nicas_randomization
  bump_nicas_subsamp_h
  bump_nicas_subsamp_hv
  bump_nicas_subsamp_hvh
  bump_nicas_subsamp_vh
  bump_obsop_adjusted
  bump_obsop_local
  bump_obsop_random)

list(APPEND saber_test_2
  bump_read-write_cmat
  bump_read-write_mom
  bump_read-write_nicas
  bump_read-write_obsop
  bump_read-write_sampling
  bump_read-write_vbal
)

# Make directories
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testdata)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref)

# Data/reference directories
foreach(test ${saber_test})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testdata/${test})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
         ${CMAKE_CURRENT_SOURCE_DIR}/testref/${test}
         ${CMAKE_CURRENT_BINARY_DIR}/testref/${test})
endforeach(test)
foreach(test ${saber_test_2})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testdata/${test}_1)
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
         ${CMAKE_CURRENT_SOURCE_DIR}/testref/${test}_1
         ${CMAKE_CURRENT_BINARY_DIR}/testref/${test}_1)
endforeach(test)

# Inter-tests links
foreach(test ${saber_test_2})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
         ${CMAKE_CURRENT_BINARY_DIR}/testdata/${test}_1
         ${CMAKE_CURRENT_BINARY_DIR}/testdata/${test}_2)
endforeach(test)
foreach(test ${saber_test_2})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_BINARY_DIR}/testref/${test}_1
        ${CMAKE_CURRENT_BINARY_DIR}/testref/${test}_2)
endforeach(test)

# Log directories
foreach(test ${saber_test})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput/${test})
endforeach(test)
foreach(test ${saber_test_2})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput/${test}_1)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput/${test}_2)
endforeach(test)


# Monocore tests: setup / run / compare
set(mpi 1)
set(omp 1)
foreach(test ${saber_test})
    ecbuild_add_test(TARGET test_${test}_${mpi}-${omp}_setup
                     TYPE SCRIPT
                     COMMAND ${CMAKE_BINARY_DIR}/bin/setup.sh
                     ARGS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${test} ${mpi} ${omp})

    ecbuild_add_test(TARGET test_${test}_${mpi}-${omp}_run
                     MPI    ${mpi}
                     OMP    ${omp}
                     COMMAND ${CMAKE_BINARY_DIR}/bin/bump.x
                     ARGS testinput/${test}_${mpi}-${omp}.nam testoutput/${test}
                     DEPENDS bump.x
                     TEST_DEPENDS test_${test}_${mpi}-${omp}_setup)

    ecbuild_add_test(TARGET test_${test}_${mpi}-${omp}_compare
                     TYPE SCRIPT
                     COMMAND ${CMAKE_BINARY_DIR}/bin/comparesaber.sh
                     ARGS ${test} ${mpi} ${omp}
                     TEST_DEPENDS test_${test}_${mpi}-${omp}_run)
endforeach()
foreach(test ${saber_test_2})
    foreach(suffix 1 2)
        ecbuild_add_test(TARGET test_${test}_${suffix}_${mpi}-${omp}_setup
                         TYPE SCRIPT
                         COMMAND ${CMAKE_BINARY_DIR}/bin/setup.sh
                         ARGS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${test}_${suffix} ${mpi} ${omp})
    endforeach()

    ecbuild_add_test(TARGET test_${test}_1_${mpi}-${omp}_run
                     MPI    ${mpi}
                     OMP    ${omp}
                     COMMAND ${CMAKE_BINARY_DIR}/bin/bump.x
                     ARGS testinput/${test}_1_${mpi}-${omp}.nam testoutput/${test}_1
                     DEPENDS bump.x
                     TEST_DEPENDS test_${test}_1_${mpi}-${omp}_setup)

    ecbuild_add_test(TARGET test_${test}_2_${mpi}-${omp}_run
                     MPI    ${mpi}
                     OMP    ${omp}
                     COMMAND ${CMAKE_BINARY_DIR}/bin/bump.x
                     ARGS testinput/${test}_2_${mpi}-${omp}.nam testoutput/${test}_2
                     DEPENDS bump.x
                     TEST_DEPENDS test_${test}_2_${mpi}-${omp}_setup test_${test}_1_${mpi}-${omp}_run)

    foreach(suffix 1 2)
        ecbuild_add_test(TARGET test_${test}_${suffix}_${mpi}-${omp}_compare
                         TYPE SCRIPT
                         COMMAND ${CMAKE_BINARY_DIR}/bin/comparesaber.sh
                         ARGS ${test}_${suffix} ${mpi} ${omp}
                         TEST_DEPENDS test_${test}_${suffix}_${mpi}-${omp}_run)
    endforeach()
endforeach()

# Dualcore test: MPI
set(SABER_MPI_TEST 1)
if(DEFINED ENV{SABER_MPI_TEST})
    set(SABER_MPI_TEST "$ENV{SABER_MPI_TEST}")
endif()
if(${SABER_MPI_TEST} MATCHES 1)
    set(mpi 2)
    set(omp 1)
    foreach(test ${saber_test})
        ecbuild_add_test(TARGET test_${test}_${mpi}-${omp}_setup
                         TYPE SCRIPT
                         COMMAND ${CMAKE_BINARY_DIR}/bin/setup.sh
                         ARGS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${test} ${mpi} ${omp}
                         TEST_DEPENDS test_${test}_1-1_setup)

        ecbuild_add_test(TARGET test_${test}_${mpi}-${omp}_run
                         MPI    ${mpi}
                         OMP    ${omp}
                         COMMAND ${CMAKE_BINARY_DIR}/bin/bump.x
                         ARGS testinput/${test}_${mpi}-${omp}.nam testoutput/${test}
                         DEPENDS bump.x
                         TEST_DEPENDS test_${test}_${mpi}-${omp}_setup)

        ecbuild_add_test(TARGET test_${test}_${mpi}-${omp}_compare
                         TYPE SCRIPT
                         COMMAND ${CMAKE_BINARY_DIR}/bin/comparesaber.sh
                         ARGS ${test} ${mpi} ${omp}
                         TEST_DEPENDS test_${test}_${mpi}-${omp}_run)
    endforeach()
    foreach(test ${saber_test_2})
        foreach(suffix 1 2)
            ecbuild_add_test(TARGET test_${test}_${suffix}_${mpi}-${omp}_setup
                             TYPE SCRIPT
                             COMMAND ${CMAKE_BINARY_DIR}/bin/setup.sh
                             ARGS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${test}_${suffix} ${mpi} ${omp})
        endforeach()

        ecbuild_add_test(TARGET test_${test}_1_${mpi}-${omp}_run
                         MPI    ${mpi}
                         OMP    ${omp}
                         COMMAND ${CMAKE_BINARY_DIR}/bin/bump.x
                         ARGS testinput/${test}_1_${mpi}-${omp}.nam testoutput/${test}_1
                         DEPENDS bump.x
                         TEST_DEPENDS test_${test}_1_${mpi}-${omp}_setup)

        ecbuild_add_test(TARGET test_${test}_2_${mpi}-${omp}_run
                         MPI    ${mpi}
                         OMP    ${omp}
                         COMMAND ${CMAKE_BINARY_DIR}/bin/bump.x
                         ARGS testinput/${test}_2_${mpi}-${omp}.nam testoutput/${test}_2
                         DEPENDS bump.x
                         TEST_DEPENDS test_${test}_2_${mpi}-${omp}_setup test_${test}_1_${mpi}-${omp}_run)

        foreach(suffix 1 2)
            ecbuild_add_test(TARGET test_${test}_${suffix}_${mpi}-${omp}_compare
                             TYPE SCRIPT
                             COMMAND ${CMAKE_BINARY_DIR}/bin/comparesaber.sh
                             ARGS ${test}_${suffix} ${mpi} ${omp}
                             TEST_DEPENDS test_${test}_${suffix}_${mpi}-${omp}_run)
        endforeach()
    endforeach()
endif()

# Dualcore test: OpenMP
set(SABER_OMP_TEST 1)
if(DEFINED ENV{SABER_OMP_TEST})
    set(SABER_OMP_TEST "$ENV{SABER_OMP_TEST}")
endif()
if(${SABER_OMP_TEST} MATCHES 1)
    set(mpi 1)
    set(omp 2)
    foreach(test ${saber_test})
        ecbuild_add_test(TARGET test_${test}_${mpi}-${omp}_setup
                         TYPE SCRIPT
                         COMMAND ${CMAKE_BINARY_DIR}/bin/setup.sh
                         ARGS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${test} ${mpi} ${omp}
                         TEST_DEPENDS test_${test}_1-1_setup)

        ecbuild_add_test(TARGET test_${test}_${mpi}-${omp}_run
                         MPI    ${mpi}
                         OMP    ${omp}
                         COMMAND ${CMAKE_BINARY_DIR}/bin/bump.x
                         ARGS testinput/${test}_${mpi}-${omp}.nam testoutput/${test}
                         DEPENDS bump.x
                         TEST_DEPENDS test_${test}_${mpi}-${omp}_setup)

        ecbuild_add_test(TARGET test_${test}_${mpi}-${omp}_compare
                         TYPE SCRIPT
                         COMMAND ${CMAKE_BINARY_DIR}/bin/comparesaber.sh
                         ARGS ${test} ${mpi} ${omp}
                         TEST_DEPENDS test_${test}_${mpi}-${omp}_run)
    endforeach()
    foreach(test ${saber_test_2})
        foreach(suffix 1 2)
            ecbuild_add_test(TARGET test_${test}_${suffix}_${mpi}-${omp}_setup
                             TYPE SCRIPT
                             COMMAND ${CMAKE_BINARY_DIR}/bin/setup.sh
                             ARGS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${test}_${suffix} ${mpi} ${omp})
        endforeach()

        ecbuild_add_test(TARGET test_${test}_1_${mpi}-${omp}_run
                         MPI    ${mpi}
                         OMP    ${omp}
                         COMMAND ${CMAKE_BINARY_DIR}/bin/bump.x
                         ARGS testinput/${test}_1_${mpi}-${omp}.nam testoutput/${test}_1
                         DEPENDS bump.x
                         TEST_DEPENDS test_${test}_1_${mpi}-${omp}_setup)

        ecbuild_add_test(TARGET test_${test}_2_${mpi}-${omp}_run
                         MPI    ${mpi}
                         OMP    ${omp}
                         COMMAND ${CMAKE_BINARY_DIR}/bin/bump.x
                         ARGS testinput/${test}_2_${mpi}-${omp}.nam testoutput/${test}_2
                         DEPENDS bump.x
                         TEST_DEPENDS test_${test}_2_${mpi}-${omp}_setup test_${test}_1_${mpi}-${omp}_run)

        foreach(suffix 1 2)
            ecbuild_add_test(TARGET test_${test}_${suffix}_${mpi}-${omp}_compare
                             TYPE SCRIPT
                             COMMAND ${CMAKE_BINARY_DIR}/bin/comparesaber.sh
                             ARGS ${test}_${suffix} ${mpi} ${omp}
                             TEST_DEPENDS test_${test}_${suffix}_${mpi}-${omp}_run)
        endforeach()
    endforeach()
endif()

# Specific compare
ecbuild_add_test(TARGET test_bump_nicas_mpicom_lsqrt_a-b_dirac
                 TYPE SCRIPT
                 COMMAND ${CMAKE_BINARY_DIR}/bin/compare_specific.sh
                 ARGS bump_nicas_mpicom_lsqrt_a bump_nicas_mpicom_lsqrt_b dirac
                 TEST_DEPENDS test_bump_nicas_mpicom_lsqrt_a_1-1_run
                              test_bump_nicas_mpicom_lsqrt_b_1-1_run)

ecbuild_add_test(TARGET test_bump_nicas_mpicom_lsqrt_a-c_dirac
                 TYPE SCRIPT
                 COMMAND ${CMAKE_BINARY_DIR}/bin/compare_specific.sh
                 ARGS bump_nicas_mpicom_lsqrt_a bump_nicas_mpicom_lsqrt_c dirac
                 TEST_DEPENDS test_bump_nicas_mpicom_lsqrt_a_1-1_run
                              test_bump_nicas_mpicom_lsqrt_c_1-1_run)
