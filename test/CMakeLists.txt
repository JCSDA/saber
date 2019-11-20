# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Test selection: TIER and CODE
set(SABER_TEST_TIER 1)
if(DEFINED ENV{SABER_TEST_TIER})
    set(SABER_TEST_TIER "$ENV{SABER_TEST_TIER}")
endif()
set(SABER_TEST_MODEL 0)
if(DEFINED ENV{SABER_TEST_MODEL})
    set(SABER_TEST_MODEL "$ENV{SABER_TEST_MODEL}")
endif()
set(SABER_TEST_MODEL_DIR "")
if(DEFINED ENV{SABER_TEST_MODEL_DIR})
    set(SABER_TEST_MODEL_DIR "$ENV{SABER_TEST_MODEL_DIR}")
endif()

# Summary
message(STATUS "SABER_TEST_TIER:      " ${SABER_TEST_TIER})
message(STATUS "SABER_TEST_MODEL:     " ${SABER_TEST_MODEL})
if(${SABER_MODEL_DIR} MATCHES 1)
    message(STATUS "-> SABER_TEST_MODEL_DIR: " ${SABER_TEST_MODEL_DIR})
endif()

# TIER 1

# Monocore tests, set 1
file(STRINGS testlist/saber_test_1.txt saber_test_tmp )
list(APPEND saber_test ${saber_test_tmp})
list(APPEND saber_data saber_data.tar.gz)
list(APPEND saber_ref saber_ref_1.tar.gz)
file(STRINGS testlist/saber_ref_to_data.txt saber_ref_to_data_tmp )
list(APPEND saber_ref_to_data ${saber_ref_to_data_tmp})

# TIER > 1
if(${SABER_TEST_TIER} GREATER 1)
    # Monocore tests, set 2
    file(STRINGS testlist/saber_test_2.txt saber_test_tmp )
    list(APPEND saber_test ${saber_test_tmp})
    list(APPEND saber_ref saber_ref_2.tar.gz)
    file(STRINGS testlist/saber_ref_to_data_mpi.txt saber_ref_to_data_tmp )
    list(APPEND saber_ref_to_data ${saber_ref_to_data_tmp})

    # Dualcore tests MPI, set 1 and 2
    file(STRINGS testlist/saber_test_1.txt saber_test_tmp )
    list(APPEND saber_test_mpi ${saber_test_tmp})
    file(STRINGS testlist/saber_test_2.txt saber_test_tmp )
    list(APPEND saber_test_mpi ${saber_test_tmp})
    list(APPEND saber_data saber_data_mpi.tar.gz)
    list(APPEND saber_ref saber_ref_mpi_1.tar.gz)
    list(APPEND saber_ref saber_ref_mpi_2.tar.gz)

    # Dualcore tests OpenMP, set 1 and 2
    file(STRINGS testlist/saber_test_1.txt saber_test_tmp )
    list(APPEND saber_test_omp ${saber_test_tmp})
    file(STRINGS testlist/saber_test_2.txt saber_test_tmp )
    list(APPEND saber_test_omp ${saber_test_tmp})
    list(APPEND saber_data saber_data_omp.tar.gz)
endif()

# TIER > 2
if(${SABER_TEST_TIER} GREATER 2)
    # Monocore tests, set 2
    file(STRINGS testlist/saber_test_3.txt saber_test_tmp )
    list(APPEND saber_test ${saber_test_tmp})
    list(APPEND saber_ref saber_ref_3.tar.gz)

    # Dualcore tests MPI, set 3
    file(STRINGS testlist/saber_test_3.txt saber_test_tmp )
    list(APPEND saber_test_mpi ${saber_test_tmp})
    list(APPEND saber_ref saber_ref_mpi_3.tar.gz)

    # Dualcore tests OpenMP, set 3
    file(STRINGS testlist/saber_test_3.txt saber_test_tmp )
    list(APPEND saber_test_omp ${saber_test_tmp})

    # Quadcore tests MPI, specific set
    file(STRINGS testlist/saber_test_quad.txt saber_test_quad_tmp )
    list(APPEND saber_test_quad ${saber_test_quad_tmp})
    list(APPEND saber_ref saber_ref_quad.tar.gz)
endif()

# TIER > 3
if(${SABER_TEST_TIER} GREATER 3)
    # Monocore tests using valgrind, set 1, 2 and 3
    file(STRINGS testlist/saber_test_1.txt saber_test_tmp )
    list(APPEND saber_test_valgrind ${saber_test_tmp})
    file(STRINGS testlist/saber_test_2.txt saber_test_tmp )
    list(APPEND saber_test_valgrind ${saber_test_tmp})
    file(STRINGS testlist/saber_test_3.txt saber_test_tmp )
    list(APPEND saber_test_valgrind ${saber_test_tmp})
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

# Set download URL
set(ECBUILD_DOWNLOAD_BASE_URL https://jedi-test-files.s3.amazonaws.com/saber)

# Download and extract data files
ecbuild_get_test_multidata(TARGET  get_saber_data
                           NAMES   ${saber_data}
                           DIRNAME testdata
                           DIRSAVE testdata
                           EXTRACT)

# Download and extract reference files
ecbuild_get_test_multidata(TARGET  get_saber_ref
                           NAMES   ${saber_ref}
                           DIRNAME testref
                           DIRSAVE testref
                           EXTRACT)

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
                     TEST_DEPENDS test_${test}_${mpi}-${omp}_run
                                  get_saber_ref)
endforeach()

# Dualcore tests MPI
set(mpi 2)
set(omp 1)
foreach(test ${saber_test_mpi})
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
                     TEST_DEPENDS test_${test}_${mpi}-${omp}_run
                                  get_saber_ref)
endforeach()

# Dualcore tests OpenMP
set(mpi 1)
set(omp 2)
foreach(test ${saber_test_omp})
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
                     TEST_DEPENDS test_${test}_${mpi}-${omp}_run
                                  get_saber_ref)
endforeach()

# Quadcore tests MPI
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
                     TEST_DEPENDS test_${test}_${mpi}-${omp}_run
                                  get_saber_ref)
endforeach()

# Monocore tests using valgrind
set(mpi 1)
set(omp 1)
foreach(test ${saber_test_valgrind})
    execute_process(COMMAND     sed "-e s/_MPI_/${mpi}/g;s/_OMP_/${omp}/g"
                    INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                    OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}_valgrind.yaml)

    ecbuild_add_test(TARGET       test_${test}_valgrind_run
                     MPI          ${mpi}
                     OMP          ${omp}
                     COMMAND      ${CMAKE_BINARY_DIR}/bin/valgrind.sh
                     ARGS         ${CMAKE_BINARY_DIR}/bin/bump.x testinput/${test}_valgrind.yaml testoutput
                     DEPENDS      bump.x
                     TEST_DEPENDS get_saber_data)
endforeach()

# Specific tests to compare results from TIER > 1
if(${SABER_TEST_TIER} GREATER 1)
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

# Model tests
if(${SABER_TEST_MODEL})
    set(mpi 4)
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
