# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

if( oops_FOUND )
    add_subdirectory( quench )
endif()

# Use find_branch function in jedi-cmake
include( git_functions )

# Default SABER_TEST_TIER
set( SABER_TEST_TIER 1 )

# Override SABER_TEST_TIER using environment variable
if( DEFINED ENV{SABER_TEST_TIER} )
    set( SABER_TEST_TIER $ENV{SABER_TEST_TIER} )
endif()

# Default test selection variables (TIER-dependent)
set( SABER_TEST_MPI 1 )
if( SABER_TEST_TIER EQUAL 1 )
    set( SABER_TEST_OMP 0 )
    set( SABER_TEST_MULTI 0 )
elseif( SABER_TEST_TIER GREATER 1 )
    set( SABER_TEST_OMP 1 )
    set( SABER_TEST_MULTI 1 )
else()
    message( FATAL_ERROR "SABER_TEST_TIER should be 1 or 2, not ${SABER_TEST_TIER}" )
endif()

# Default test selection variables (TIER-independent)
set( SABER_TEST_VALGRIND 0 )
set( SABER_TEST_BUMP_QUENCH 0 )
set( SABER_TEST_GSI_GFS 0 )
set( SABER_TEST_GSI_GEOS 0 )
set( SABER_TEST_SPECTRALB 0 )
if( oops_FOUND )
    set( SABER_TEST_BUMP_QUENCH 1 )
    if( gsibec_FOUND )
        set( SABER_TEST_GSI_GEOS 1 )
    endif()
    if( gsibec_FOUND AND sp_FOUND )
        set( SABER_TEST_GSI_GFS 1 )
    endif()
    if( ENABLE_SPECTRALB )
        set( SABER_TEST_SPECTRALB 1 )
    endif()
endif()

# Override test selection variables using environment variables
if( DEFINED ENV{SABER_TEST_MPI} )
    set( SABER_TEST_MPI $ENV{SABER_TEST_MPI} )
endif()
if( DEFINED ENV{SABER_TEST_OMP} )
    set( SABER_TEST_OMP $ENV{SABER_TEST_OMP} )
endif()
if( DEFINED ENV{SABER_TEST_MULTI} )
    set( SABER_TEST_MULTI $ENV{SABER_TEST_MULTI} )
endif()
if( DEFINED ENV{SABER_TEST_VALGRIND} )
    set( SABER_TEST_VALGRIND $ENV{SABER_TEST_VALGRIND} )
endif()
if( oops_FOUND )
    if( DEFINED ENV{SABER_TEST_BUMP_QUENCH} )
        set( SABER_TEST_BUMP_QUENCH $ENV{SABER_TEST_BUMP_QUENCH} )
    endif()
    if( gsibec_FOUND )
        if( DEFINED ENV{SABER_TEST_GSI_GEOS} )
            set( SABER_TEST_GSI_GEOS $ENV{SABER_TEST_GSI_GEOS} )
        endif()
    endif()
    if( gsibec_FOUND AND sp_FOUND )
        if( DEFINED ENV{SABER_TEST_GSI_GFS} )
            set( SABER_TEST_GSI_GFS $ENV{SABER_TEST_GSI_GFS} )
        endif()
    endif()
    if( ENABLE_SPECTRALB )
        if( DEFINED ENV{SABER_TEST_SPECTRALB} )
            set( SABER_TEST_SPECTRALB $ENV{SABER_TEST_SPECTRALB} )
        endif()
    endif()
endif()

# Summary
message( STATUS "SABER variables:" )
message( STATUS "  - SABER_TEST_TIER:          ${SABER_TEST_TIER}" )
message( STATUS "  - SABER_TEST_MPI:           ${SABER_TEST_MPI}" )
message( STATUS "  - SABER_TEST_OMP:           ${SABER_TEST_OMP}" )
message( STATUS "  - SABER_TEST_MULTI:         ${SABER_TEST_MULTI}" )
message( STATUS "  - SABER_TEST_VALGRIND:      ${SABER_TEST_VALGRIND}" )
message( STATUS "  - SABER_TEST_BUMP_QUENCH:   ${SABER_TEST_BUMP_QUENCH}" )
message( STATUS "  - SABER_TEST_GSI_GEOS:      ${SABER_TEST_GSI_GEOS}" )
message( STATUS "  - SABER_TEST_GSI_GFS:       ${SABER_TEST_GSI_GFS}" )
message( STATUS "  - SABER_TEST_SPECTRALB:     ${SABER_TEST_SPECTRALB}" )

# Tests
message( STATUS "SABER tests activated:" )

# Input data
file( STRINGS testlist/saber_data.txt saber_data_tmp )
list( APPEND saber_data ${saber_data_tmp} )
list( APPEND saber_data_tar saber_data.tar.gz)

# TIER 1
message( STATUS "  - TIER 1 base" )
file( STRINGS testlist/saber_test_tier1.txt saber_test_tmp )
list( APPEND saber_test ${saber_test_tmp} )
list( APPEND saber_test_full ${saber_test_tmp} )
file( STRINGS testlist/saber_test_tier1-cgal.txt saber_test_cgal_tmp )
list( APPEND saber_test_full ${saber_test_cgal_tmp} )
if( atlas_HAVE_TESSELATION )
    message( STATUS "  - TIER 1 CGAL-specific" )
    list( APPEND saber_test ${saber_test_cgal_tmp} )
endif()
file( STRINGS testlist/saber_test_tier1-multi.txt saber_test_multi_tmp )
list( APPEND saber_test_full ${saber_test_multi_tmp} )
if( SABER_TEST_MULTI )
    message( STATUS "  - TIER 1 multicores" )
    list( APPEND saber_test_multi ${saber_test_multi_tmp} )
endif()
file( STRINGS testlist/saber_test_tier1-bump-quench.txt saber_test_bump_quench_tmp )
list( APPEND saber_test_full ${saber_test_bump_quench_tmp} )
if( SABER_TEST_BUMP_QUENCH )
    message( STATUS "  - TIER 1 BUMP-QUENCH-specific" )
    list( APPEND saber_test_bump_quench ${saber_test_bump_quench_tmp} )
endif()
file( STRINGS testlist/saber_test_tier1-gsi-geos.txt saber_test_gsi_geos_tmp )
list( APPEND saber_test_full ${saber_test_gsi_geos_tmp} )
if( SABER_TEST_GSI_GEOS )
    message( STATUS "  - TIER 1 GSI-GEOS-specific" )
    list( APPEND saber_test_gsi_geos ${saber_test_gsi_geos_tmp} )
endif()
file( STRINGS testlist/saber_test_tier1-gsi-gfs.txt saber_test_gsi_gfs_tmp )
list( APPEND saber_test_full ${saber_test_gsi_gfs_tmp} )
if( SABER_TEST_GSI_GFS )
    message( STATUS "  - TIER 1 GSI-GFS-specific" )
    list( APPEND saber_test_gsi_gfs ${saber_test_gsi_gfs_tmp} )
endif()

file( STRINGS testlist/saber_test_tier1-spectralb.txt saber_test_spectralb_tmp )
list( APPEND saber_test_full ${saber_test_spectralb_tmp} )
if( SABER_TEST_SPECTRALB )
    message( STATUS "  - TIER 1 SPECTRALB-specific" )
    list( APPEND saber_test_spectralb ${saber_test_spectralb_tmp} )
endif()
file( STRINGS testlist/saber_ref_tier1.txt saber_ref_tmp )
list( APPEND saber_ref ${saber_ref_tmp} )
list( APPEND saber_ref_tar saber_ref_tier1.tar.gz)

# TIER > 1
if( SABER_TEST_TIER GREATER 1 )
    message( STATUS "  - TIER 2 base" )
    file( STRINGS testlist/saber_test_tier2.txt saber_test_tmp )
    list( APPEND saber_test_full ${saber_test_tmp} )
    list( APPEND saber_test ${saber_test_tmp} )
    file( STRINGS testlist/saber_ref_tier2.txt saber_ref_tmp )
    list( APPEND saber_ref ${saber_ref_tmp} )
    list( APPEND saber_ref_tar saber_ref_tier2.tar.gz )
endif()

# Set download URL

# Define branch name
if( DEFINED ENV{SABER_TESTFILES_BRANCH} )
    # Branch specified by environment
    set( GIT_BRANCH_SABER $ENV{SABER_TESTFILES_BRANCH} )
    message( STATUS "Branch name provided by user: ${GIT_BRANCH_SABER}" )
else()
    # Find branch name
    find_branch_name( REPO_DIR_NAME saber )

    if( DEFINED GIT_TAG_FUNC )
        # Tag is checked out => download data from DASH
        set( GIT_BRANCH_SABER ${GIT_TAG_FUNC} )
        message( STATUS "SABER is in tag: ${GIT_BRANCH_SABER}" )
    elseif( DEFINED GIT_BRANCH_FUNC )
        # Branch is checked out => data already downloaded in saber-data
        set( GIT_BRANCH_SABER ${GIT_BRANCH_FUNC} )
        message( STATUS "SABER is in branch: ${GIT_BRANCH_SABER}" )
    endif()
endif()

# Define where the data are already present or need to be downloaded (and how)
if( DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES} )
    # Data should be present in the path defined by LOCAL_PATH_JEDI_TESTFILES
    if( EXISTS $ENV{LOCAL_PATH_JEDI_TESTFILES}/saber/${GIT_BRANCH_SABER} )
        set( SABER_TESTFILES_PATH $ENV{LOCAL_PATH_JEDI_TESTFILES}/saber/${GIT_BRANCH_SABER} )
    else()
        set( SABER_TESTFILES_PATH $ENV{LOCAL_PATH_JEDI_TESTFILES}/saber/develop )
    endif()
    message( STATUS "Files loaded from: " ${SABER_TESTFILES_PATH} )
elseif( DEFINED GIT_BRANCH_FUNC )
    # Data loaded from saber-data repo
    set( SABER_TESTFILES_PATH ${CMAKE_SOURCE_DIR}/saber-data )
    message( STATUS "Files loaded from: " ${SABER_TESTFILES_PATH} )
else( DEFINED GIT_TAG_FUNC )
    # ECBUILD_DOWNLOAD_BASE_URL env var can be used to force test files and
    # crtm coef to be downloaded from other databases such as S3 instead of DASH
    # example ECBUILD_DOWNLOAD_BASE_URL=https://jedi-test-files.s3.amazonaws.com

    if( DEFINED ENV{ECBUILD_DOWNLOAD_BASE_URL} )
        set( ECBUILD_DOWNLOAD_BASE_URL "$ENV{ECBUILD_DOWNLOAD_BASE_URL}/saber/${GIT_TAG_FUNC}" )
        set( checksum "0")
    else()
        set( ECBUILD_DOWNLOAD_BASE_URL https://gdex.ucar.edu/dataset/jedi-skylab/file )
        set( checksum "0")
    endif()
    message( STATUS "Files downloaded from: " ${ECBUILD_DOWNLOAD_BASE_URL} )

    # Set SABER_TESTFILES_PATH
    set( SABER_TESTFILES_PATH ${CMAKE_SOURCE_DIR}/test-data-release/saber/${GIT_TAG_FUNC} )
    message( STATUS "Files downloaded into: " ${SABER_TESTFILES_PATH} )

    list( APPEND SABER_DATA_DOWNLOADER_ARGS
          ${ECBUILD_DOWNLOAD_BASE_URL}
          ${CMAKE_SOURCE_DIR}/test-data-release
          "saber_testinput_tier_1_1.3.0.tar.gz"
          ${checksum} )

    # Create test-data-release in source directory
    file( MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/test-data-release )

    # Create download script
    set ( FILENAME saber_data_downloader.py )
    set ( SOURCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME} )
    set ( DEST_FILE ${CMAKE_BINARY_DIR}/bin/${FILENAME} )

    # Set configure file
    if( EXISTS "${SOURCE_FILE}.in" )
      configure_file( ${SOURCE_FILE}.in ${DEST_FILE} @ONLY )
    else()
      configure_file( ${SOURCE_FILE}    ${DEST_FILE} @ONLY )
    endif()

    # Add download script
    add_custom_target( bin_saber_test_download_script ALL
                       COMMAND chmod +x ${DEST_FILE}
                       DEPENDS ${DEST_FILE} )

    # Add test for downloading data
    ecbuild_add_test( TARGET    get_saber_test_data
                      TYPE      SCRIPT
                      COMMAND   ${CMAKE_BINARY_DIR}/bin/saber_data_downloader.py
                      ARGS  ${SABER_DATA_DOWNLOADER_ARGS} )
endif()

# Setup SABER directories and links
message( STATUS "Setup SABER directories and links" )
file( WRITE ${CMAKE_BINARY_DIR}/bin/saber_testdir )
foreach( test ${saber_test_full} )
    file( APPEND ${CMAKE_BINARY_DIR}/bin/saber_testdir "${test}\n" )
endforeach()
file( WRITE ${CMAKE_BINARY_DIR}/bin/saber_testdata )
foreach( file ${saber_data} )
    file( APPEND ${CMAKE_BINARY_DIR}/bin/saber_testdata "${file}\n" )
endforeach()
file( WRITE ${CMAKE_BINARY_DIR}/bin/saber_testref )
foreach( file ${saber_ref} )
    file( APPEND ${CMAKE_BINARY_DIR}/bin/saber_testref "${file}\n" )
endforeach()
execute_process( COMMAND bash ${CMAKE_BINARY_DIR}/bin/saber_setup.sh
                         ${SABER_TESTFILES_PATH}
                         ${CMAKE_CURRENT_BINARY_DIR}
                         ${CMAKE_BINARY_DIR} )

# Executables
ecbuild_add_executable( TARGET  saber_bump.x
                        SOURCES mains/bump_main.cc
                                mains/bump_main.fypp
                        LIBS    saber )

# Process fypp files
fckit_target_preprocess_fypp( saber_bump.x )

if( vader_FOUND )
    set( vader_LIBRARIES "vader" )
else()
    set( vader_LIBRARIES "" )
endif()
if( gsibec_FOUND )
    message( STATUS "SABER TEST LIB gsibec was found" )
    set( gsibec_LIBRARIES "gsibec" )
else()
    message( STATUS "SABER TEST LIB gsibec was NOT found" )
    set( gsibec_LIBRARIES "" )
endif()

if ( SABER_TEST_BUMP_QUENCH OR SABER_TEST_GSI_GEOS OR SABER_TEST_GSI_GFS OR SABER_TEST_SPECTRALB )
    ecbuild_add_executable( TARGET  saber_quench_randomization.x
                            SOURCES mains/quenchRandomization.cc
                            LIBS    quench
                                    ${vader_LIBRARIES}
                                    ${gsibec_LIBRARIES}
                                    saber )

    ecbuild_add_executable( TARGET  saber_quench_error_covariance_training.x
                            SOURCES mains/quenchErrorCovarianceTraining.cc
                            LIBS    quench
                                    ${vader_LIBRARIES}
                                    ${gsibec_LIBRARIES}
                                    saber )

    ecbuild_add_executable( TARGET  saber_quench_saber_block_test.x
                            SOURCES mains/quenchSaberBlockTest.cc
                            LIBS    quench
                                    ${vader_LIBRARIES}
                                    ${gsibec_LIBRARIES}
                                    saber )

    ecbuild_add_executable( TARGET  saber_quench_dirac.x
                            SOURCES mains/quenchDirac.cc
                            LIBS    quench
                                    ${vader_LIBRARIES}
                                    ${gsibec_LIBRARIES}
                                    saber )

    ecbuild_add_executable( TARGET  saber_quench_convertstate.x
                            SOURCES mains/quenchConvertState.cc
                            LIBS    quench
                                    ${vader_LIBRARIES}
                                    saber )
endif()

# Mono-core tests
set( mpi 1 )
set( omp 1 )
foreach( test ${saber_test} )
    execute_process( COMMAND     sed "-e s/_MPI_/${mpi}/g;s/_OMP_/${omp}/g"
                     INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                     OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}_${mpi}-${omp}.yaml )

    ecbuild_add_test( TARGET       saber_test_${test}_${mpi}-${omp}_run
                      MPI          ${mpi}
                      OMP          ${omp}
                      COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_bump.x
                      ARGS         testinput/${test}_${mpi}-${omp}.yaml testoutput
                      DEPENDS      saber_bump.x
                      TEST_DEPENDS get_saber_data )

    ecbuild_add_test( TARGET       saber_test_${test}_${mpi}-${omp}_compare
                      TYPE SCRIPT
                      COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                      ARGS         ${test} ${mpi} ${omp}
                      TEST_DEPENDS get_saber_ref
                                   saber_test_${test}_${mpi}-${omp}_run )

    if( SABER_TEST_VALGRIND )
        ecbuild_add_test( TARGET saber_test_${test}_${mpi}-${omp}_valgrind
                          MPI    ${mpi}
                          OMP    ${omp}
                          COMMAND ${CMAKE_BINARY_DIR}/bin/saber_valgrind.sh
                          ARGS ${CMAKE_BINARY_DIR}/bin/saber_bump.x testinput/${test}_${mpi}-${omp}.yaml testoutput
                          DEPENDS saber_bump.x
                          TEST_DEPENDS get_saber_data )
    endif()
endforeach()

# Dual-core tests (MPI)
if( SABER_TEST_MPI )
    set( mpi 2 )
    set( omp 1 )
    foreach( test ${saber_test} )
        execute_process( COMMAND     sed "-e s/_MPI_/${mpi}/g;s/_OMP_/${omp}/g"
                         INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                         OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}_${mpi}-${omp}.yaml )

        ecbuild_add_test( TARGET       saber_test_${test}_${mpi}-${omp}_run
                          MPI          ${mpi}
                          OMP          ${omp}
                          COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_bump.x
                          ARGS         testinput/${test}_${mpi}-${omp}.yaml testoutput
                          DEPENDS      saber_bump.x
                          TEST_DEPENDS get_saber_data )

        ecbuild_add_test( TARGET       saber_test_${test}_${mpi}-${omp}_compare
                          TYPE SCRIPT
                          COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                          ARGS         ${test} ${mpi} ${omp}
                          TEST_DEPENDS get_saber_ref
                                       saber_test_${test}_${mpi}-${omp}_run )
    endforeach()
endif()

# Dual-core tests (OpenMP)
if( SABER_TEST_OMP )
    set( mpi 1 )
    set( omp 2 )
    foreach( test ${saber_test} )
        execute_process( COMMAND     sed "-e s/_MPI_/${mpi}/g;s/_OMP_/${omp}/g"
                         INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                         OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}_${mpi}-${omp}.yaml )

        ecbuild_add_test( TARGET       saber_test_${test}_${mpi}-${omp}_run
                          MPI          ${mpi}
                          OMP          ${omp}
                          COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_bump.x
                          ARGS         testinput/${test}_${mpi}-${omp}.yaml testoutput
                          DEPENDS      saber_bump.x
                          TEST_DEPENDS get_saber_data )

        ecbuild_add_test( TARGET       saber_test_${test}_${mpi}-${omp}_compare
                          TYPE SCRIPT
                          COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                          ARGS         ${test} ${mpi} ${omp}
                          TEST_DEPENDS get_saber_ref
                                       saber_test_${test}_${mpi}-${omp}_run )
    endforeach()
endif()

# Multi-core tests (MPI)
if( SABER_TEST_MULTI )
    set( omp 1 )
    foreach( mpi RANGE 4 12 2 )
        foreach( test ${saber_test_multi} )
            execute_process( COMMAND     sed "-e s/_MPI_/${mpi}/g;s/_OMP_/${omp}/g"
                             INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                             OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}_${mpi}-${omp}.yaml )

            ecbuild_add_test( TARGET       saber_test_${test}_${mpi}-${omp}_run
                              MPI          ${mpi}
                              OMP          ${omp}
                              COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_bump.x
                              ARGS         testinput/${test}_${mpi}-${omp}.yaml testoutput
                              DEPENDS      saber_bump.x
                              TEST_DEPENDS get_saber_data )

            ecbuild_add_test( TARGET       saber_test_${test}_${mpi}-${omp}_compare
                              TYPE SCRIPT
                              COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                              ARGS         ${test} ${mpi} ${omp}
                              TEST_DEPENDS get_saber_ref
                                           saber_test_${test}_${mpi}-${omp}_run )
        endforeach()
    endforeach()
endif()

# Post-comparisons (between tests)

# Compare read and write tests
list( APPEND read-write_test mom )
list( APPEND read-write_test nicas_global )
list( APPEND read-write_test nicas_local )
list( APPEND read-write_test sampling_global )
list( APPEND read-write_test sampling_local )
list( APPEND read-write_test var )
list( APPEND read-write_test vbal_cov_global )
list( APPEND read-write_test vbal_cov_local )
list( APPEND read-write_test vbal_global )
list( APPEND read-write_test vbal_local )
foreach( test ${read-write_test} )
    # Mono-core tests
    set( mpi 1 )
    set( omp 1 )
    ecbuild_add_test( TARGET       saber_test_bump_read-write_${test}_${mpi}-${omp}_post
                      TYPE SCRIPT
                      COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                      ARGS         bump_read_${test} bump_write_${test} ${mpi}-${omp}
                      TEST_DEPENDS saber_test_bump_read_${test}_${mpi}-${omp}_run
                                   saber_test_bump_write_${test}_${mpi}-${omp}_run )

    # Dual-core tests (MPI)
    if( SABER_TEST_MPI )
        set( mpi 2 )
        set( omp 1 )
        ecbuild_add_test( TARGET       saber_test_bump_read-write_${test}_${mpi}-${omp}_post
                          TYPE SCRIPT
                          COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                          ARGS         bump_read_${test} bump_write_${test} ${mpi}-${omp}
                          TEST_DEPENDS saber_test_bump_read_${test}_${mpi}-${omp}_run
                                       saber_test_bump_write_${test}_${mpi}-${omp}_run )
    endif()

    # Dual-core tests (OpenMP)
    if( SABER_TEST_OMP )
        set( mpi 1 )
        set( omp 2 )
        ecbuild_add_test( TARGET       saber_test_bump_read-write_${test}_${mpi}-${omp}_post
                          TYPE SCRIPT
                          COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                          ARGS         bump_read_${test} bump_write_${test} ${mpi}-${omp}
                          TEST_DEPENDS saber_test_bump_read_${test}_${mpi}-${omp}_run
                                       saber_test_bump_write_${test}_${mpi}-${omp}_run )
    endif()

    # Multi-core tests (MPI)
    if( SABER_TEST_MULTI )
        set( omp 1 )
        foreach( mpi RANGE 4 12 2 )
            ecbuild_add_test( TARGET       saber_test_bump_read-write_${test}_${mpi}-${omp}_post
                              TYPE SCRIPT
                              COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                              ARGS         bump_read_${test} bump_write_${test} ${mpi}-${omp}
                              TEST_DEPENDS saber_test_bump_read_${test}_${mpi}-${omp}_run
                                           saber_test_bump_write_${test}_${mpi}-${omp}_run )
        endforeach()
    endif()
endforeach()

if( SABER_TEST_TIER GREATER 1 )
    # Compare ensemble-based and covariance-based vertical balance
    ecbuild_add_test( TARGET       saber_test_bump_vbal_ens-cov_post
                      TYPE SCRIPT
                      COMMAND      ${CMAKE_BINARY_DIR}/bin/saber_compare.sh
                      ARGS         bump_vbal_no_mask bump_vbal_no_mask_cov 1-1
                      TEST_DEPENDS saber_test_bump_vbal_no_mask_1-1_run
                                   saber_test_bump_vbal_no_mask_cov_1-1_run )
endif()

# OOPS-based tests

# BUMP-QUENCH tests
if( SABER_TEST_BUMP_QUENCH )
    # 2 MPI tasks and 2 OpenMP threads
    set( mpi 2 )
    set( omp 2 )

    # Link to yaml files
    foreach( test ${saber_test_bump_quench} )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                         ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                         ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}.yaml )
    endforeach()

    # Randomization tests
    foreach( test ${saber_test_bump_quench} )
        string( FIND ${test} "quench_randomization" start_index )
        if( start_index MATCHES 0 )
            ecbuild_add_test( TARGET saber_test_${test}
                              MPI ${mpi}
                              OMP ${omp}
                              COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_randomization.x
                              ARGS testinput/${test}.yaml
                              DEPENDS saber_quench_randomization.x )
        endif()
    endforeach()

    # Error covariance training tests
    foreach( test ${saber_test_bump_quench} )
        string( FIND ${test} "quench_error_covariance_training" start_index )
        if( start_index MATCHES 0 )
            ecbuild_add_test( TARGET saber_test_${test}
                              MPI ${mpi}
                              OMP ${omp}
                              COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_training.x
                              ARGS testinput/${test}.yaml
                              DEPENDS saber_quench_error_covariance_training.x )
        endif()
    endforeach()

    # Adjoint/inverse tests
    foreach( test ${saber_test_bump_quench} )
        string( FIND ${test} "quench_saber_block_test" start_index )
        if( start_index MATCHES 0 )
            ecbuild_add_test( TARGET saber_test_${test}
                              MPI ${mpi}
                              OMP ${omp}
                              COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_saber_block_test.x
                              ARGS testinput/${test}.yaml
                              DEPENDS saber_quench_saber_block_test.x )
        endif()
    endforeach()

    # Dirac tests
    foreach( test ${saber_test_bump_quench} )
        string( FIND ${test} "quench_dirac" start_index )
        if( start_index MATCHES 0 )
            ecbuild_add_test( TARGET saber_test_${test}
                              MPI ${mpi}
                              OMP ${omp}
                              COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_dirac.x
                              ARGS testinput/${test}.yaml
                              DEPENDS saber_quench_dirac.x )
        endif()
    endforeach()
endif()

if( SABER_TEST_GSI_GEOS )
    # 1 MPI task and 2 OpenMP threads
    set( mpi 1 )
    set( omp 2 )

    # Link to yaml files
    foreach( test ${saber_test_gsi_geos} )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                         ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                         ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}.yaml )
    endforeach()

    # Link to nml files
    foreach( test ${saber_test_gsi_geos} )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                         ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.nml
                         ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}.nml )
    endforeach()

    # Dirac tests
    foreach( test ${saber_test_gsi_geos} )
        string( FIND ${test} "quench_dirac" start_index )
        if( start_index MATCHES 0 )
            ecbuild_add_test( TARGET saber_test_${test}
                              MPI ${mpi}
                              OMP ${omp}
                              COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_dirac.x
                              ARGS testinput/${test}.yaml
                              DEPENDS saber_quench_dirac.x )
        endif()
    endforeach()
endif()

if ( SABER_TEST_GSI_GFS )
    # 1 MPI task and 2 OpenMP threads
    set( mpi 1 )
    set( omp 2 )

    # Link to yaml files
    foreach( test ${saber_test_gsi_gfs} )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                         ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                         ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}.yaml )
     endforeach()

    # Link to nml files
    foreach( test ${saber_test_gsi_gfs} )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                         ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.nml
                         ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}.nml )
    endforeach()

    # Dirac tests
    foreach( test ${saber_test_gsi_gfs} )
        string( FIND ${test} "quench_dirac" start_index )
        if( start_index MATCHES 0 )
            ecbuild_add_test( TARGET saber_test_${test}
                              MPI ${mpi}
                              OMP ${omp}
                              COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_dirac.x
                              ARGS testinput/${test}.yaml
                              DEPENDS saber_quench_dirac.x )
        endif()
    endforeach()
endif()

if( SABER_TEST_SPECTRALB )
    # 1 OpenMP thread
    set( omp 1 )

    # Link to yaml files
    foreach( test ${saber_test_spectralb} )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                         ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                         ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}.yaml )
    endforeach()

    # Adjoint/inverse tests
    foreach( test ${saber_test_spectralb} )
        string( FIND ${test} "quench_saber_block_test" start_index )
        if( start_index MATCHES 0 )
            if ( test MATCHES quench_saber_block_test_spectralb_from_CS )
                set( mpi 1 )
            elseif(  test MATCHES quench_saber_block_test_spectralb_from_L15  )
                set( mpi 1 )
            else()
                set( mpi 4 )
            endif ()
            ecbuild_add_test( TARGET saber_test_${test}
                              MPI ${mpi}
                              OMP ${omp}
                              COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_saber_block_test.x
                              ARGS testinput/${test}.yaml
                              DEPENDS saber_quench_saber_block_test.x )
        endif()
    endforeach()

    # Dirac tests
    foreach( test ${saber_test_spectralb} )
        string( FIND ${test} "quench_dirac" start_index )
        if( start_index MATCHES 0 )
            set( mpi 1 )
            ecbuild_add_test( TARGET saber_test_${test}
                              MPI ${mpi}
                              OMP ${omp}
                              COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_dirac.x
                              ARGS testinput/${test}.yaml
                              DEPENDS saber_quench_dirac.x )
        endif()
    endforeach()
endif()
