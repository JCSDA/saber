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
    set(SABER_TEST_QUAD 0)
    set(SABER_TEST_OMP 1)
elseif(${SABER_TEST_TIER} MATCHES 3)
    set(SABER_TEST_MPI 1)
    set(SABER_TEST_QUAD 1)
    set(SABER_TEST_OMP 1)
else()
    message(FATAL_ERROR "SABER_TEST_TIER should be 1, 2 or 3, not " ${SABER_TEST_TIER})
endif()

# Default test selection variables (TIER-independent)
set(SABER_TEST_VALGRIND 0)
set(SABER_TEST_MODEL 0)
set(SABER_TEST_MODEL_DIR "")

# Override test selection variables using environment variables
if(DEFINED ENV{SABER_TEST_MPI})
    set(SABER_TEST_MPI "$ENV{SABER_TEST_MPI}")
endif()
if(DEFINED ENV{SABER_TEST_QUAD})
    set(SABER_TEST_QUAD "$ENV{SABER_TEST_QUAD}")
endif()
if(DEFINED ENV{SABER_TEST_OMP})
    set(SABER_TEST_OMP "$ENV{SABER_TEST_OMP}")
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

# Summary
message(STATUS "SABER_TEST_TIER:      " ${SABER_TEST_TIER})
message(STATUS "SABER_TEST_MPI:       " ${SABER_TEST_MPI})
message(STATUS "SABER_TEST_OMP:       " ${SABER_TEST_OMP})
message(STATUS "SABER_TEST_QUAD:      " ${SABER_TEST_QUAD})
message(STATUS "SABER_TEST_VALGRIND:  " ${SABER_TEST_VALGRIND})
message(STATUS "SABER_TEST_MODEL:     " ${SABER_TEST_MODEL})
if(${SABER_TEST_MODEL} MATCHES 1)
    message(STATUS "  SABER_TEST_MODEL_DIR: " ${SABER_TEST_MODEL_DIR})
endif()

# TIER 1
file(STRINGS testlist/saber_test_1.txt saber_test_tmp )
list(APPEND saber_test ${saber_test_tmp})
file(STRINGS testlist/saber_data.txt saber_data_tmp )
list(APPEND saber_data ${saber_data_tmp})
list(APPEND saber_data_tar saber_data.tar.gz)
file(STRINGS testlist/saber_ref_1.txt saber_ref_tmp )
list(APPEND saber_ref ${saber_ref_tmp})
list(APPEND saber_ref_tar saber_ref_1.tar.gz)
file(STRINGS testlist/saber_ref_to_data.txt saber_ref_to_data_tmp )
list(APPEND saber_ref_to_data ${saber_ref_to_data_tmp})
if(${SABER_TEST_MPI})
    file(STRINGS testlist/saber_data_mpi.txt saber_data_tmp )
    list(APPEND saber_data ${saber_data_tmp})
    list(APPEND saber_data_tar saber_data_mpi.tar.gz)
    file(STRINGS testlist/saber_ref_mpi_1.txt saber_ref_tmp )
    list(APPEND saber_ref ${saber_ref_tmp})
    list(APPEND saber_ref_tar saber_ref_mpi_1.tar.gz)
    file(STRINGS testlist/saber_ref_to_data_mpi.txt saber_ref_to_data_tmp )
    list(APPEND saber_ref_to_data ${saber_ref_to_data_tmp})
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
if(${SABER_TEST_QUAD} MATCHES 1)
    file(STRINGS testlist/saber_test_quad.txt saber_test_quad_tmp )
    list(APPEND saber_test_quad ${saber_test_quad_tmp})
    file(STRINGS testlist/saber_ref_quad.txt saber_ref_tmp )
    list(APPEND saber_ref ${saber_ref_tmp})
    list(APPEND saber_ref_tar saber_ref_quad.tar.gz)
endif()

# Model tests
if(${SABER_TEST_MODEL} MATCHES 1)
    file(STRINGS testlist/saber_test_model.txt saber_test_model_tmp )
    list(APPEND saber_test_model ${saber_test_model_tmp})
    file(STRINGS testlist/saber_data_model.txt saber_data_tmp )
    list(APPEND saber_data ${saber_data_tmp})
    list(APPEND saber_data_tar saber_data_model.tar.gz)
endif()

# Make directories
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testdata)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref)

# Set download URL
set(ECBUILD_DOWNLOAD_BASE_URL https://jedi-test-files.s3.amazonaws.com/saber)

# Download and extract data files
ecbuild_get_test_multidata(TARGET  get_saber_data
                           NAMES   ${saber_data_tar}
                           DIRNAME "."
                           DIRSAVE testdata
                           EXTRACT)

# Link data files into build directory
foreach(file ${saber_data})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                    ${PROJECT_SOURCE_DIR}/test_data_download/testdata/${file}
                    ${CMAKE_CURRENT_BINARY_DIR}/testdata/${file})
endforeach()
foreach(file ${saber_ref_to_data})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                    ${PROJECT_SOURCE_DIR}/test_data_download/testref/${file}
                    ${CMAKE_CURRENT_BINARY_DIR}/testdata/${file})
endforeach()

# Download and extract reference files
ecbuild_get_test_multidata(TARGET  get_saber_ref
                           NAMES   ${saber_ref_tar}
                           DIRNAME "."
                           DIRSAVE testref
                           EXTRACT)

# Link reference files into build directory
foreach(file ${saber_ref})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                    ${PROJECT_SOURCE_DIR}/test_data_download/testref/${file}
                    ${CMAKE_CURRENT_BINARY_DIR}/testref/${file})
endforeach()

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
if(${SABER_TEST_MPI} MATCHES 1)
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
if(${SABER_TEST_OMP} MATCHES 1)
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
if(${SABER_TEST_QUAD} MATCHES 1)
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
