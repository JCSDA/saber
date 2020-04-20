# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Default SABER_TEST_TIER
set(SABER_TEST_TIER 1)

# Override SABER_TEST_TIER using environment variable
if(DEFINED ENV{SABER_TEST_TIER})
    set(SABER_TEST_TIER "$ENV{SABER_TEST_TIER}")
endif()

# Default test selection variables (TIER-dependent)
if(${SABER_TEST_TIER} MATCHES 1)
    set(SABER_TEST_MPI 1)
    set(SABER_TEST_OMP 0)
    set(SABER_TEST_QUAD 0)
elseif(${SABER_TEST_TIER} MATCHES 2)
    set(SABER_TEST_MPI 1)
    set(SABER_TEST_OMP 1)
    set(SABER_TEST_QUAD 0)
elseif(${SABER_TEST_TIER} MATCHES 3)
    set(SABER_TEST_MPI 1)
    set(SABER_TEST_OMP 1)
    set(SABER_TEST_QUAD 1)
else()
    message(FATAL_ERROR "SABER_TEST_TIER should be 1, 2 or 3, not " ${SABER_TEST_TIER})
endif()

# Default test selection variables (TIER-independent)
set(SABER_TEST_VALGRIND 0)
set(SABER_TEST_MODEL 0)
set(SABER_TEST_MODEL_DIR "")
if (${OOPS_FOUND})
   set(SABER_TEST_OOPS 1)
else()
   set(SABER_TEST_OOPS 0)
endif()

# Override test selection variables using environment variables
if(DEFINED ENV{SABER_TEST_MPI})
    set(SABER_TEST_MPI "$ENV{SABER_TEST_MPI}")
endif()
if(DEFINED ENV{SABER_TEST_OMP})
    set(SABER_TEST_OMP "$ENV{SABER_TEST_OMP}")
endif()
if(DEFINED ENV{SABER_TEST_QUAD})
    set(SABER_TEST_QUAD "$ENV{SABER_TEST_QUAD}")
endif()
if(DEFINED ENV{SABER_TEST_VALGRIND})
    set(SABER_TEST_VALGRIND "$ENV{SABER_TEST_VALGRIND}")
endif()
if(DEFINED ENV{SABER_TEST_MODEL})
    set(SABER_TEST_MODEL "$ENV{SABER_TEST_MODEL}")
endif()
if(DEFINED ENV{SABER_TEST_MODEL_DIR})
    set(SABER_TEST_MODEL_DIR "$ENV{SABER_TEST_MODEL_DIR}")
endif()
if(DEFINED ENV{SABER_TEST_OOPS})
    set(SABER_TEST_OOPS "$ENV{SABER_TEST_OOPS}")
endif()

# Summary
message(STATUS "SABER_TEST_TIER:      " ${SABER_TEST_TIER})
message(STATUS "SABER_TEST_MPI:       " ${SABER_TEST_MPI})
message(STATUS "SABER_TEST_OMP:       " ${SABER_TEST_OMP})
message(STATUS "SABER_TEST_QUAD:      " ${SABER_TEST_QUAD})
message(STATUS "SABER_TEST_VALGRIND:  " ${SABER_TEST_VALGRIND})
message(STATUS "SABER_TEST_MODEL:     " ${SABER_TEST_MODEL})
if(${SABER_TEST_MODEL})
    message(STATUS "  SABER_TEST_MODEL_DIR: " ${SABER_TEST_MODEL_DIR})
endif()
message(STATUS "SABER_TEST_OOPS:      " ${SABER_TEST_OOPS})

# TIER 1
file(STRINGS testlist/saber_test_1.txt saber_test_tmp )
list(APPEND saber_test ${saber_test_tmp})
file(STRINGS testlist/saber_data.txt saber_data_tmp )
list(APPEND saber_data ${saber_data_tmp})
list(APPEND saber_data_tar saber_data.tar.gz)
file(STRINGS testlist/saber_ref_1.txt saber_ref_tmp )
list(APPEND saber_ref ${saber_ref_tmp})
list(APPEND saber_ref_tar saber_ref_1.tar.gz)
if(${SABER_TEST_MPI})
    file(STRINGS testlist/saber_data_mpi.txt saber_data_tmp )
    list(APPEND saber_data ${saber_data_tmp})
    list(APPEND saber_data_tar saber_data_mpi.tar.gz)
    file(STRINGS testlist/saber_ref_mpi_1.txt saber_ref_tmp )
    list(APPEND saber_ref ${saber_ref_tmp})
    list(APPEND saber_ref_tar saber_ref_mpi_1.tar.gz)
endif()
if(${SABER_TEST_OMP})
    file(STRINGS testlist/saber_data_omp.txt saber_data_tmp )
    list(APPEND saber_data ${saber_data_tmp})
    list(APPEND saber_data_tar saber_data_omp.tar.gz)
endif()

# TIER > 1
if(${SABER_TEST_TIER} GREATER 1)
    file(STRINGS testlist/saber_test_2.txt saber_test_tmp )
    list(APPEND saber_test ${saber_test_tmp})
    file(STRINGS testlist/saber_ref_2.txt saber_ref_tmp )
    list(APPEND saber_ref ${saber_ref_tmp})
    list(APPEND saber_ref_tar saber_ref_2.tar.gz)
    if(${SABER_TEST_MPI})
        file(STRINGS testlist/saber_ref_mpi_2.txt saber_ref_tmp )
        list(APPEND saber_ref ${saber_ref_tmp})
        list(APPEND saber_ref_tar saber_ref_mpi_2.tar.gz)
    endif()
endif()

# TIER > 2
if(${SABER_TEST_TIER} GREATER 2)
    file(STRINGS testlist/saber_test_3.txt saber_test_tmp )
    list(APPEND saber_test ${saber_test_tmp})
    file(STRINGS testlist/saber_ref_3.txt saber_ref_tmp )
    list(APPEND saber_ref ${saber_ref_tmp})
    list(APPEND saber_ref_tar saber_ref_3.tar.gz)
    if(${SABER_TEST_MPI})
        file(STRINGS testlist/saber_ref_mpi_3.txt saber_ref_tmp )
        list(APPEND saber_ref ${saber_ref_tmp})
        list(APPEND saber_ref_tar saber_ref_mpi_3.tar.gz)
    endif()
endif()

# Quad-core tests
if(${SABER_TEST_QUAD})
    file(STRINGS testlist/saber_test_quad.txt saber_test_quad_tmp )
    list(APPEND saber_test_quad ${saber_test_quad_tmp})
    file(STRINGS testlist/saber_ref_quad.txt saber_ref_tmp )
    list(APPEND saber_ref ${saber_ref_tmp})
    list(APPEND saber_ref_tar saber_ref_quad.tar.gz)
endif()

# Model tests
if(${SABER_TEST_MODEL})
    file(STRINGS testlist/saber_test_model.txt saber_test_model_tmp )
    list(APPEND saber_test_model ${saber_test_model_tmp})
endif()

# OOPS tests
if(${SABER_TEST_OOPS})
    file(STRINGS testlist/saber_test_oops.txt saber_test_oops_tmp )
    list(APPEND saber_test_oops ${saber_test_oops_tmp})
    file(STRINGS testlist/saber_data_oops.txt saber_data_tmp )
    list(APPEND saber_data ${saber_data_tmp})
    list(APPEND saber_data_tar saber_data_oops.tar.gz)
    file(STRINGS testlist/saber_ref_oops.txt saber_ref_tmp )
    list(APPEND saber_ref ${saber_ref_tmp})
    list(APPEND saber_ref_tar saber_ref_oops.tar.gz)
endif()

# Set download URL
set(ECBUILD_DOWNLOAD_BASE_URL https://jedi-test-files.s3.amazonaws.com/saber)

# Define local path for test files
if(DEFINED ENV{LOCAL_PATH_TESTFILES_SABER})
    set(LOCAL_PATH_TESTFILES_SABER "$ENV{LOCAL_PATH_TESTFILES_SABER}")
endif()

# If local path to testfiles is defined don't download
if( NOT DEFINED LOCAL_PATH_TESTFILES_SABER )
    # Set TESTFILE_DIR_SABER
    set(TESTFILE_DIR_SABER ${CMAKE_BINARY_DIR}/test_data/saber CACHE PATH "data dir for test data")

    # Get the current git branch
    execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    message(STATUS "Files downloaded into: " ${TESTFILE_DIR_SABER} "/" ${GIT_BRANCH})
    message(STATUS "Download sources:")

    # Make specific data directory
    file(MAKE_DIRECTORY ${TESTFILE_DIR_SABER}/${GIT_BRANCH})

    # Check whether the URL exists or not
    ecbuild_check_multiurl(NAMES   ${saber_data_tar}
                           DIRHOST ${GIT_BRANCH}
                           RESULT  SPECIFIC_TEST_FILES)

    # Set distant directory
    if(${SPECIFIC_TEST_FILES} MATCHES 0)
        # Download and extract new test files (distant directory = git branch)
        set(DIRNAME ${GIT_BRANCH})
    else()
        # Download and extract develop test files (distant directory = develop)
        set(DIRNAME "develop")
    endif()
    message(STATUS "  Test data: " ${ECBUILD_DOWNLOAD_BASE_URL} "/" ${DIRNAME})

    # Download and extract data files
    ecbuild_get_test_multidata(TARGET   get_saber_data
                               NAMES    ${saber_data_tar}
                               DIRNAME  ${DIRNAME}
                               DIRLOCAL ${TESTFILE_DIR_SABER}/${GIT_BRANCH}/testdata
                               EXTRACT)

    # Check whether the URL exists or not
    ecbuild_check_multiurl(NAMES   ${saber_ref_tar}
                           DIRHOST ${GIT_BRANCH}
                           RESULT  SPECIFIC_TEST_FILES)

    # Set distant directory
    if(${SPECIFIC_TEST_FILES} MATCHES 0)
        # Download and extract new test files (distant directory = git branch)
        set(DIRNAME ${GIT_BRANCH})
    else()
        # Download and extract develop test files (distant directory = develop)
        set(DIRNAME "develop")
    endif()
    message(STATUS "  Reference: " ${ECBUILD_DOWNLOAD_BASE_URL}/${DIRNAME})

    # Download and extract reference files
    ecbuild_get_test_multidata(TARGET   get_saber_ref
                               NAMES    ${saber_ref_tar}
                               DIRNAME  ${DIRNAME}
                               DIRLOCAL ${TESTFILE_DIR_SABER}/${GIT_BRANCH}/testref
                               EXTRACT)
else()
    set(TESTFILE_DIR_SABER ${LOCAL_PATH_TESTFILES_SABER})
    message(STATUS "Files loaded from: " ${TESTFILE_DIR_SABER})
endif()

# Make directories
file(MAKE_DIRECTORY ${TESTFILE_DIR_SABER}/${GIT_BRANCH}/testdata)
file(MAKE_DIRECTORY ${TESTFILE_DIR_SABER}/${GIT_BRANCH}/testref)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testdata)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref)
foreach(test ${saber_test})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testdata/${test})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput/${test})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref/${test})
endforeach()
foreach(test ${saber_test_model})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testdata/${test})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput/${test})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref/${test})
endforeach()
foreach(test ${saber_test_oops})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testdata/${test})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput/${test})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref/${test})
endforeach()

# Link data/reference files
foreach(file ${saber_data})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                    ${TESTFILE_DIR_SABER}/${GIT_BRANCH}/testdata/${file}
                    ${CMAKE_CURRENT_BINARY_DIR}/testdata/${file})
endforeach()
foreach(file ${saber_ref})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                    ${TESTFILE_DIR_SABER}/${GIT_BRANCH}/testref/${file}
                    ${CMAKE_CURRENT_BINARY_DIR}/testref/${file})
endforeach()

if(NOT ${METIS_FOUND})
    # Remove METIS tests
    file(STRINGS testlist/saber_test_metis.txt saber_test_tmp )
    list(REMOVE_ITEM saber_test ${saber_test_tmp})
    if(${SABER_TEST_QUAD})
        file(STRINGS testlist/saber_test_metis.txt saber_test_quad_tmp )
        list(REMOVE_ITEM saber_test_quad ${saber_test_quad_tmp})
    endif()
endif()

# Mono-core tests
set(mpi 1)
set(omp 1)
foreach(test ${saber_test})
    execute_process(COMMAND     sed "-e s/_MPI_/${mpi}/g;s/_OMP_/${omp}/g"
                    INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                    OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}_${mpi}-${omp}.yaml)

    ecbuild_add_test(TARGET       test_${test}_${mpi}-${omp}_run
                     MPI          ${mpi}
                     OMP          ${omp}
                     COMMAND      ${CMAKE_BINARY_DIR}/bin/bump.x
                     ARGS         testinput/${test}_${mpi}-${omp}.yaml testoutput
                     DEPENDS      bump.x
                     TEST_DEPENDS get_saber_data)
 
    ecbuild_add_test(TARGET       test_${test}_${mpi}-${omp}_compare
                     TYPE SCRIPT
                     COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                     ARGS         ${test} ${mpi} ${omp}
                     TEST_DEPENDS get_saber_ref
                                  test_${test}_${mpi}-${omp}_run)

    if(${SABER_TEST_VALGRIND})
        ecbuild_add_test(TARGET test_${test}_${mpi}-${omp}_valgrind
                         MPI    ${mpi}
                         OMP    ${omp}
                         COMMAND ${CMAKE_BINARY_DIR}/bin/valgrind.sh
                         ARGS ${CMAKE_BINARY_DIR}/bin/bump.x testinput/${test}_${mpi}-${omp}.yaml testoutput
                         DEPENDS bump.x
                         TEST_DEPENDS get_saber_data)
    endif()
endforeach()

# Dual-core tests (MPI)
if(${SABER_TEST_MPI})
    set(mpi 2)
    set(omp 1)
    foreach(test ${saber_test})
        execute_process(COMMAND     sed "-e s/_MPI_/${mpi}/g;s/_OMP_/${omp}/g"
                        INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                        OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}_${mpi}-${omp}.yaml)

        ecbuild_add_test(TARGET       test_${test}_${mpi}-${omp}_run
                         MPI          ${mpi}
                         OMP          ${omp}
                         COMMAND      ${CMAKE_BINARY_DIR}/bin/bump.x
                         ARGS         testinput/${test}_${mpi}-${omp}.yaml testoutput
                         DEPENDS      bump.x
                         TEST_DEPENDS get_saber_data)
   
        ecbuild_add_test(TARGET       test_${test}_${mpi}-${omp}_compare
                         TYPE SCRIPT
                         COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                         ARGS         ${test} ${mpi} ${omp}
                         TEST_DEPENDS get_saber_ref
                                      test_${test}_${mpi}-${omp}_run)
    endforeach()
endif()

# Dual-core tests (OpenMP)
if(${SABER_TEST_OMP})
    set(mpi 1)
    set(omp 2)
    foreach(test ${saber_test})
        execute_process(COMMAND     sed "-e s/_MPI_/${mpi}/g;s/_OMP_/${omp}/g"
                        INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                        OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}_${mpi}-${omp}.yaml)

        ecbuild_add_test(TARGET       test_${test}_${mpi}-${omp}_run
                         MPI          ${mpi}
                         OMP          ${omp}
                         COMMAND      ${CMAKE_BINARY_DIR}/bin/bump.x
                         ARGS         testinput/${test}_${mpi}-${omp}.yaml testoutput
                         DEPENDS      bump.x
                         TEST_DEPENDS get_saber_data)
    
        ecbuild_add_test(TARGET       test_${test}_${mpi}-${omp}_compare
                         TYPE SCRIPT
                         COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                         ARGS         ${test} ${mpi} ${omp}
                         TEST_DEPENDS get_saber_ref
                                      test_${test}_${mpi}-${omp}_run)
    endforeach()
endif()

# Quad-core tests (MPI)
if(${SABER_TEST_QUAD})
    set(mpi 4)
    set(omp 1)
    foreach(test ${saber_test_quad})
        execute_process(COMMAND     sed "-e s/_MPI_/${mpi}/g;s/_OMP_/${omp}/g"
                        INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                        OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}_${mpi}-${omp}.yaml)

        ecbuild_add_test(TARGET       test_${test}_${mpi}-${omp}_run
                         MPI          ${mpi}
                         OMP          ${omp}
                         COMMAND      ${CMAKE_BINARY_DIR}/bin/bump.x
                         ARGS         testinput/${test}_${mpi}-${omp}.yaml testoutput
                         DEPENDS      bump.x
                         TEST_DEPENDS get_saber_data)

        ecbuild_add_test(TARGET       test_${test}_${mpi}-${omp}_compare
                         TYPE SCRIPT
                         COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                         ARGS         ${test} ${mpi} ${omp}
                         TEST_DEPENDS get_saber_ref
                                      test_${test}_${mpi}-${omp}_run)
    endforeach()
endif()

# Specific compare (TIER > 1)
if(${SABER_TEST_TIER} GREATER 1)
    ecbuild_add_test(TARGET       test_bump_nicas_mpicom_lsqrt_a-b_dirac
                     TYPE SCRIPT
                     COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                     ARGS         bump_nicas_mpicom_lsqrt_a bump_nicas_mpicom_lsqrt_b dirac
                     TEST_DEPENDS test_bump_nicas_mpicom_lsqrt_a_1-1_run
                                  test_bump_nicas_mpicom_lsqrt_b_1-1_run)

    ecbuild_add_test(TARGET       test_bump_nicas_mpicom_lsqrt_a-c_dirac
                     TYPE SCRIPT
                     COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                     ARGS         bump_nicas_mpicom_lsqrt_a bump_nicas_mpicom_lsqrt_c dirac
                     TEST_DEPENDS test_bump_nicas_mpicom_lsqrt_a_1-1_run
                                  test_bump_nicas_mpicom_lsqrt_c_1-1_run)
endif()

# Model tests
if(${SABER_TEST_MODEL})
    set(mpi 6)
    set(omp 1)
    foreach(test ${saber_test_model})
        execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                        ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                        ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}.yaml)

        ecbuild_add_test(TARGET       test_${test}_link
                         COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_links.ksh
                         ARGS         ${SABER_TEST_MODEL_DIR} ${CMAKE_CURRENT_BINARY_DIR}/testdata ${test})

        ecbuild_add_test(TARGET       test_${test}_run
                         MPI          ${mpi}
                         OMP          ${omp}
                         COMMAND      ${CMAKE_BINARY_DIR}/bin/bump.x
                         ARGS         testinput/${test}.yaml testoutput
                         DEPENDS      bump.x
                         TEST_DEPENDS test_${test}_link)
    endforeach()
endif()

# OOPS tests
if(${SABER_TEST_OOPS})
    #--------------------------------------------------------------------
    # Link to yaml files
    #--------------------------------------------------------------------

    foreach(test ${saber_test_oops})
        execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                        ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                        ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}.yaml)
    endforeach()

    if (${ENABLE_OOPS_TOYMODELS})
        #--------------------------------------------------------------------
        # Parameters tests
        #--------------------------------------------------------------------

        ecbuild_add_test( TARGET test_qg_parameters_bump_cov_run
                          MPI    1
                          OMP    2
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_qg_estimate_parameters.x
                          ARGS testinput/qg_parameters_bump_cov.yaml testoutput/qg_parameters_bump_cov/test.log.out
                          DEPENDS saber_qg_estimate_parameters.x
                          TEST_DEPENDS get_saber_data )

        #--------------------------------------------------------------------

        ecbuild_add_test( TARGET test_qg_parameters_bump_hyb_run
                          MPI    1
                          OMP    2
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_qg_estimate_parameters.x
                          ARGS testinput/qg_parameters_bump_hyb.yaml testoutput/qg_parameters_bump_hyb/test.log.out
                          DEPENDS saber_qg_estimate_parameters.x
                          TEST_DEPENDS get_saber_data )

        #--------------------------------------------------------------------

        ecbuild_add_test( TARGET test_qg_parameters_bump_lct_run
                          MPI    1
                          OMP    2
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_qg_estimate_parameters.x
                          ARGS testinput/qg_parameters_bump_lct.yaml testoutput/qg_parameters_bump_lct/test.log.out
                          DEPENDS saber_qg_estimate_parameters.x
                          TEST_DEPENDS get_saber_data )

        #--------------------------------------------------------------------

        ecbuild_add_test( TARGET test_qg_parameters_bump_loc_3d_run
                          MPI    1
                          OMP    2
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_qg_estimate_parameters.x
                          ARGS testinput/qg_parameters_bump_loc_3d.yaml testoutput/qg_parameters_bump_loc_3d/test.log.out
                          DEPENDS saber_qg_estimate_parameters.x
                          TEST_DEPENDS get_saber_data )

        #--------------------------------------------------------------------

        ecbuild_add_test( TARGET test_qg_parameters_bump_loc_4d_run
                          MPI    1
                          OMP    2
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_qg_estimate_parameters.x
                          ARGS testinput/qg_parameters_bump_loc_4d.yaml testoutput/qg_parameters_bump_loc_4d/test.log.out
                          DEPENDS saber_qg_estimate_parameters.x
                          TEST_DEPENDS get_saber_data )

        #--------------------------------------------------------------------
        # BUMP links
        #--------------------------------------------------------------------

        ecbuild_add_test( TARGET test_qg_bump_links
                          TYPE SCRIPT
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_oops_links.sh
                          TEST_DEPENDS test_qg_parameters_bump_cov_run
                                       test_qg_parameters_bump_hyb_run
                                       test_qg_parameters_bump_lct_run
                                       test_qg_parameters_bump_loc_3d_run
                                       test_qg_parameters_bump_loc_4d_run )

        #--------------------------------------------------------------------
        # BUMP dirac tests
        #--------------------------------------------------------------------

        ecbuild_add_test( TARGET test_qg_dirac_bump_cov_run
                          MPI    1
                          OMP    2
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_qg_dirac.x
                          ARGS testinput/qg_dirac_bump_cov.yaml testoutput/qg_dirac_bump_cov/test.log.out
                          DEPENDS saber_qg_dirac.x
                          TEST_DEPENDS test_qg_bump_links )

        #--------------------------------------------------------------------

        ecbuild_add_test( TARGET test_qg_dirac_bump_hyb_run
                          MPI    1
                          OMP    2
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_qg_dirac.x
                          ARGS testinput/qg_dirac_bump_hyb.yaml testoutput/qg_dirac_bump_hyb/test.log.out
                          DEPENDS saber_qg_dirac.x
                          TEST_DEPENDS test_qg_bump_links )

        #--------------------------------------------------------------------

        ecbuild_add_test( TARGET test_qg_dirac_bump_lct_run
                          MPI    1
                          OMP    2
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_qg_dirac.x
                          ARGS testinput/qg_dirac_bump_lct.yaml testoutput/qg_dirac_bump_lct/test.log.out
                          DEPENDS saber_qg_dirac.x
                          TEST_DEPENDS test_qg_bump_links )

        #--------------------------------------------------------------------

        ecbuild_add_test( TARGET test_qg_dirac_bump_loc_3d_run
                          MPI    1
                          OMP    2
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_qg_dirac.x
                          ARGS testinput/qg_dirac_bump_loc_3d.yaml testoutput/qg_dirac_bump_loc_3d/test.log.out
                          DEPENDS saber_qg_dirac.x
                          TEST_DEPENDS test_qg_bump_links )

        #--------------------------------------------------------------------

        ecbuild_add_test( TARGET test_qg_dirac_bump_loc_4d_run
                          MPI    1
                          OMP    2
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_qg_dirac.x
                          ARGS testinput/qg_dirac_bump_loc_4d.yaml testoutput/qg_dirac_bump_loc_4d/test.log.out
                          DEPENDS saber_qg_dirac.x
                          TEST_DEPENDS test_qg_bump_links )

        #--------------------------------------------------------------------
        # 3d variational tests
        #--------------------------------------------------------------------

        ecbuild_add_test( TARGET test_qg_3densvar_bump_run
                          MPI    1
                          OMP    2
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_qg_4dvar.x
                          ARGS testinput/qg_3densvar_bump.yaml testoutput/qg_3densvar_bump/test.log.out
                          DEPENDS saber_qg_4dvar.x
                          TEST_DEPENDS test_qg_bump_links )

        #--------------------------------------------------------------------

        ecbuild_add_test( TARGET test_qg_3dvar_bump_run
                          MPI    1
                          OMP    2
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_qg_4dvar.x
                          ARGS testinput/qg_3dvar_bump.yaml testoutput/qg_3dvar_bump/test.log.out
                          DEPENDS saber_qg_4dvar.x
                          TEST_DEPENDS test_qg_bump_links )

        #--------------------------------------------------------------------

        ecbuild_add_test( TARGET test_qg_3dvar_hybrid_bump_run
                          MPI    1
                          OMP    2
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_qg_4dvar.x
                          ARGS testinput/qg_3dvar_hybrid_bump.yaml testoutput/qg_3dvar_hybrid_bump/test.log.out
                          DEPENDS saber_qg_4dvar.x
                          TEST_DEPENDS test_qg_bump_links )

        #--------------------------------------------------------------------
        # 4d variational tests
        #--------------------------------------------------------------------

        ecbuild_add_test( TARGET test_qg_4densvar_bump_run
                          MPI    1
                          OMP    2
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_qg_4dvar.x
                          ARGS testinput/qg_4densvar_bump.yaml testoutput/qg_4densvar_bump/test.log.out
                          DEPENDS saber_qg_4dvar.x
                          TEST_DEPENDS test_qg_bump_links )

        #--------------------------------------------------------------------

        ecbuild_add_test( TARGET test_qg_4dvar_drplanczos_bump_run
                          MPI    1
                          OMP    2
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_qg_4dvar.x
                          ARGS testinput/qg_4dvar_drplanczos_bump.yaml testoutput/qg_4dvar_drplanczos_bump/test.log.out
                          DEPENDS saber_qg_4dvar.x
                          TEST_DEPENDS test_qg_bump_links )

        #--------------------------------------------------------------------

        ecbuild_add_test( TARGET test_qg_4dvar_drplanczos_hybrid_bump_run
                          MPI    1
                          OMP    2
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_qg_4dvar.x
                          ARGS testinput/qg_4dvar_drplanczos_hybrid_bump.yaml testoutput/qg_4dvar_drplanczos_hybrid_bump/test.log.out
                          DEPENDS saber_qg_4dvar.x
                          TEST_DEPENDS test_qg_bump_links )

        #--------------------------------------------------------------------
        # Compare tests
        #--------------------------------------------------------------------
        foreach(test ${saber_test_oops})
            ecbuild_add_test( TARGET test_${test}_compare
                              TYPE SCRIPT
                              COMMAND ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                              ARGS ${test}
                              TEST_DEPENDS get_saber_ref test_${test}_run )
        endforeach()

    endif()
endif()