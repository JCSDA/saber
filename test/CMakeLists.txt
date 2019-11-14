# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Test selection: TIER and CODE
set(SABER_TEST_TIER 1)
if(DEFINED ENV{SABER_TEST_TIER})
    set(SABER_TEST_TIER "$ENV{SABER_TEST_TIER}")
endif()
set(SABER_TEST_CODE 00000)
if(DEFINED ENV{SABER_TEST_CODE})
    set(SABER_TEST_CODE "$ENV{SABER_TEST_CODE}")
endif()
set(SABER_TEST_RUN 1)
if(DEFINED ENV{SABER_TEST_RUN})
    set(SABER_TEST_RUN "$ENV{SABER_TEST_RUN}")
endif()
set(SABER_TEST_COMPARE 1)
if(DEFINED ENV{SABER_TEST_COMPARE})
    set(SABER_TEST_COMPARE "$ENV{SABER_TEST_COMPARE}")
endif()

# Split SABER_TEST_CODE into other tests variables
set(CODE ${SABER_TEST_CODE})

# First bit: SABER_TEST_MPI
math(EXPR TMP "${CODE} / 10")
math(EXPR SABER_TEST_MPI "${CODE}-${TMP}*10")
set(CODE ${TMP})

# Second bit: SABER_TEST_MPI4
math(EXPR TMP "${CODE} / 10")
math(EXPR SABER_TEST_MPI4 "${CODE}-${TMP}*10")
set(CODE ${TMP})

# Third bit: SABER_TEST_OMP
math(EXPR TMP "${CODE} / 10")
math(EXPR SABER_TEST_OMP "${CODE}-${TMP}*10")
set(CODE ${TMP})

# Fourth bit: SABER_TEST_VALGRIND
math(EXPR TMP "${CODE} / 10")
math(EXPR SABER_TEST_VALGRIND "${CODE}-${TMP}*10")
set(CODE ${TMP})

# Fourth bit: SABER_TEST_MODEL
math(EXPR TMP "${CODE} / 10")
math(EXPR SABER_TEST_MODEL "${CODE}-${TMP}*10")
set(CODE ${TMP})

# Summary
message(STATUS "SABER_TEST_TIER: " ${SABER_TEST_TIER})
message(STATUS "SABER_TEST_CODE: " ${SABER_TEST_CODE})
message(STATUS "  SABER_TEST_MPI:       " ${SABER_TEST_MPI})
message(STATUS "  SABER_TEST_MPI4:      " ${SABER_TEST_MPI4})
message(STATUS "  SABER_TEST_OMP:       " ${SABER_TEST_OMP})
message(STATUS "  SABER_TEST_VALGRIND:  " ${SABER_TEST_VALGRIND})
message(STATUS "  SABER_TEST_MODEL:     " ${SABER_TEST_MODEL})

# TIER 1
file(STRINGS testlist/saber_test_1.txt saber_test_tmp )
list(APPEND saber_test ${saber_test_tmp})
list(APPEND saber_data saber_data.tar.gz)
list(APPEND saber_ref saber_ref_1.tar.gz)
file(STRINGS testlist/saber_ref_to_data.txt saber_ref_to_data_tmp )
list(APPEND saber_ref_to_data ${saber_ref_to_data_tmp})
if(${SABER_TEST_MPI})
    list(APPEND saber_ref saber_ref_mpi_1.tar.gz)
    file(STRINGS testlist/saber_ref_to_data_mpi.txt saber_ref_to_data_tmp )
    list(APPEND saber_ref_to_data ${saber_ref_to_data_tmp})
endif()
if(${SABER_TEST_OMP})
    file(STRINGS testlist/saber_ref_to_data_omp.txt saber_ref_to_data_tmp )
    list(APPEND saber_ref_to_data ${saber_ref_to_data_tmp})
endif()

# TIER > 1
if(${SABER_TEST_TIER} GREATER 1)
    file(STRINGS testlist/saber_test_2.txt saber_test_tmp )
    list(APPEND saber_test ${saber_test_tmp})
    list(APPEND saber_ref saber_ref_2.tar.gz)
    if(${SABER_TEST_MPI})
        list(APPEND saber_ref saber_ref_mpi_2.tar.gz)
    endif()
endif()

# TIER > 2
if(${SABER_TEST_TIER} GREATER 2)
    file(STRINGS testlist/saber_test_3.txt saber_test_tmp )
    list(APPEND saber_test ${saber_test_tmp})
    list(APPEND saber_ref saber_ref_3.tar.gz)
    if(${SABER_TEST_MPI})
        list(APPEND saber_ref saber_ref_mpi_3.tar.gz)
    endif()
endif()

# Quadcore tests (should be among saber_test)
if(${SABER_TEST_MPI4} MATCHES 1)
    file(STRINGS testlist/saber_test_quad.txt saber_test_quad_tmp )
    list(APPEND saber_test_quad ${saber_test_quad_tmp})
    list(APPEND saber_ref saber_ref_quad.tar.gz)
endif()

# Model tests
if(${SABER_TEST_MODEL} MATCHES 1)
    file(STRINGS testlist/saber_test_model.txt saber_test_model_tmp )
    list(APPEND saber_test_model ${saber_test_model_tmp})
    list(APPEND saber_data saber_data_model.tar.gz)
endif()

# Make directories
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testdata)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref)

# Download and extract data files
ecbuild_get_test_multidata(TARGET  get_saber_data
                           URL     https://jedi-test-files.s3.amazonaws.com/saber
                           NAMES   ${saber_data}
                           DIRNAME testdata
                           DIRSAVE testdata
                           EXTRACT
                           NOCHECK)

# Clean data files
ecbuild_add_test(TARGET  clean_saber_data
                         COMMAND ${CMAKE_COMMAND}
                         ARGS -E chdir ${CMAKE_CURRENT_BINARY_DIR}/testdata rm -f ${saber_data}
                         TEST_DEPENDS get_saber_data)

# Download and extract reference files
ecbuild_get_test_multidata(TARGET  get_saber_ref
                           URL     https://jedi-test-files.s3.amazonaws.com/saber
                           NAMES   ${saber_ref}
                           DIRNAME testref
                           DIRSAVE testref
                           EXTRACT
                           NOCHECK)

# Clean reference files
ecbuild_add_test(TARGET  clean_saber_ref
                         COMMAND ${CMAKE_COMMAND}
                         ARGS -E chdir ${CMAKE_CURRENT_BINARY_DIR}/testref rm -f ${saber_ref}
                         TEST_DEPENDS get_saber_ref)

# Link reference files into data directory
foreach(file ${saber_ref_to_data})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                    ${CMAKE_CURRENT_BINARY_DIR}/testref/${file}
                    ${CMAKE_CURRENT_BINARY_DIR}/testdata/${file})
endforeach()

# Monocore tests
set(mpi 1)
set(omp 1)
foreach(test ${saber_test})
    execute_process(COMMAND     sed "-e s/_MPI_/${mpi}/g;s/_OMP_/${omp}/g"
                    INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                    OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}_${mpi}-${omp}.yaml)

    if(${SABER_TEST_RUN})
        ecbuild_add_test(TARGET       test_${test}_${mpi}-${omp}_run
                         MPI          ${mpi}
                         OMP          ${omp}
                         COMMAND      ${CMAKE_BINARY_DIR}/bin/bump.x
                         ARGS         testinput/${test}_${mpi}-${omp}.yaml testoutput
                         DEPENDS      bump.x
                         TEST_DEPENDS test_${test}_${mpi}-${omp}_setup
                                      clean_saber_data)
    endif()

    if(${SABER_TEST_COMPARE})
        ecbuild_add_test(TARGET       test_${test}_${mpi}-${omp}_compare
                         TYPE SCRIPT
                         COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                         ARGS         ${test} ${mpi} ${omp}
                         TEST_DEPENDS test_${test}_${mpi}-${omp}_run
                                      clean_saber_ref)
    endif()

    if(${SABER_TEST_VALGRIND})
        ecbuild_add_test(TARGET test_${test}_${mpi}-${omp}_valgrind
                         MPI    ${mpi}
                         OMP    ${omp}
                         COMMAND ${CMAKE_BINARY_DIR}/bin/valgrind.sh
                         ARGS ${CMAKE_BINARY_DIR}/bin/bump.x testinput/${test}_${mpi}-${omp}.yaml testoutput
                         DEPENDS bump.x
                         TEST_DEPENDS test_${test}_${mpi}-${omp}_setup)
    endif()
endforeach()

# Dualcore tests (MPI)
if(${SABER_TEST_MPI} MATCHES 1)
    set(mpi 2)
    set(omp 1)
    foreach(test ${saber_test})
        execute_process(COMMAND     sed "-e s/_MPI_/${mpi}/g;s/_OMP_/${omp}/g"
                        INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                        OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}_${mpi}-${omp}.yaml)

        if(${SABER_TEST_RUN})
            ecbuild_add_test(TARGET       test_${test}_${mpi}-${omp}_run
                             MPI          ${mpi}
                             OMP          ${omp}
                             COMMAND      ${CMAKE_BINARY_DIR}/bin/bump.x
                             ARGS         testinput/${test}_${mpi}-${omp}.yaml testoutput
                             DEPENDS      bump.x
                             TEST_DEPENDS test_${test}_${mpi}-${omp}_setup
                                          clean_saber_data)
        endif()
    
        if(${SABER_TEST_COMPARE})
            ecbuild_add_test(TARGET       test_${test}_${mpi}-${omp}_compare
                             TYPE SCRIPT
                             COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                             ARGS         ${test} ${mpi} ${omp}
                             TEST_DEPENDS test_${test}_${mpi}-${omp}_run
                                          clean_saber_ref)
        endif()
    endforeach()
endif()

# Dualcore tests (OpenMP)
if(${SABER_TEST_OMP} MATCHES 1)
    set(mpi 1)
    set(omp 2)
    foreach(test ${saber_test})
        execute_process(COMMAND     sed "-e s/_MPI_/${mpi}/g;s/_OMP_/${omp}/g"
                        INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                        OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}_${mpi}-${omp}.yaml)

        if(${SABER_TEST_RUN})
            ecbuild_add_test(TARGET       test_${test}_${mpi}-${omp}_run
                             MPI          ${mpi}
                             OMP          ${omp}
                             COMMAND      ${CMAKE_BINARY_DIR}/bin/bump.x
                             ARGS         testinput/${test}_${mpi}-${omp}.yaml testoutput
                             DEPENDS      bump.x
                             TEST_DEPENDS test_${test}_${mpi}-${omp}_setup
                                          clean_saber_data)
        endif()
    
        if(${SABER_TEST_COMPARE})
            ecbuild_add_test(TARGET       test_${test}_${mpi}-${omp}_compare
                             TYPE SCRIPT
                             COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                             ARGS         ${test} ${mpi} ${omp}
                             TEST_DEPENDS test_${test}_${mpi}-${omp}_run
                                          clean_saber_ref)
        endif()
    endforeach()
endif()

# Quadcore tests (MPI)
if(${SABER_TEST_MPI4} MATCHES 1)
    set(mpi 4)
    set(omp 1)
    foreach(test ${saber_test_quad})
        execute_process(COMMAND     sed "-e s/_MPI_/${mpi}/g;s/_OMP_/${omp}/g"
                        INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                        OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}_${mpi}-${omp}.yaml)

        if(${SABER_TEST_RUN})
            ecbuild_add_test(TARGET       test_${test}_${mpi}-${omp}_run
                             MPI          ${mpi}
                             OMP          ${omp}
                             COMMAND      ${CMAKE_BINARY_DIR}/bin/bump.x
                             ARGS         testinput/${test}_${mpi}-${omp}.yaml testoutput
                             DEPENDS      bump.x
                             TEST_DEPENDS test_${test}_${mpi}-${omp}_setup
                                          clean_saber_data)
        endif()
    
        if(${SABER_TEST_COMPARE})
            ecbuild_add_test(TARGET       test_${test}_${mpi}-${omp}_compare
                             TYPE SCRIPT
                             COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                             ARGS         ${test} ${mpi} ${omp}
                             TEST_DEPENDS test_${test}_${mpi}-${omp}_run
                                          clean_saber_ref)
        endif()
    endforeach()
endif()

# Specific compare (TIER > 1)
if(${SABER_TEST_TIER} GREATER 1)
    if(${SABER_TEST_COMPARE})
        ecbuild_add_test(TARGET       test_bump_nicas_mpicom_lsqrt_a-b_dirac
                         TYPE SCRIPT
                         COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare_specific.sh
                         ARGS         bump_nicas_mpicom_lsqrt_a bump_nicas_mpicom_lsqrt_b dirac
                         TEST_DEPENDS test_bump_nicas_mpicom_lsqrt_a_1-1_run
                                      test_bump_nicas_mpicom_lsqrt_b_1-1_run)

        ecbuild_add_test(TARGET       test_bump_nicas_mpicom_lsqrt_a-c_dirac
                         TYPE SCRIPT
                         COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare_specific.sh
                         ARGS         bump_nicas_mpicom_lsqrt_a bump_nicas_mpicom_lsqrt_c dirac
                         TEST_DEPENDS test_bump_nicas_mpicom_lsqrt_a_1-1_run
                                      test_bump_nicas_mpicom_lsqrt_c_1-1_run)
    endif()
endif()

# Model tests
if(${SABER_TEST_MODEL})
    set(mpi 4)
    set(omp 1)
    foreach(test ${saber_test_model})
        execute_process(COMMAND     sed "-e s/_MPI_/${mpi}/g;s/_OMP_/${omp}/g"
                        INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                        OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}_${mpi}-${omp}.yaml)

        ecbuild_add_test(TARGET       test_${test}_${mpi}-${omp}_run
                         MPI          ${mpi}
                         OMP          ${omp}
                         COMMAND      ${CMAKE_BINARY_DIR}/bin/bump.x
                         ARGS         testinput/${test}_${mpi}-${omp}.yaml testoutput
                         DEPENDS      bump.x
                         TEST_DEPENDS test_${test}_${mpi}-${omp}_setup
                                          clean_saber_data)
    endforeach()
endif()
