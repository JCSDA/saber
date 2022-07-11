# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

if( oops_FOUND )
    # Add quench testbed
    add_subdirectory( quench )

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
    set( SABER_TEST_GSI 0 )
    if( gsibclim_FOUND)
        set( SABER_TEST_GSI 1 )
    endif()
    set( SABER_TEST_SPECTRALB 0 )
    if( atlas_TRANS_FOUND)
        set( SABER_TEST_SPECTRALB 1 )
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
    if( gsibclim_FOUND )
        if( DEFINED ENV{SABER_TEST_GSI} )
            set( SABER_TEST_GSI $ENV{SABER_TEST_GSI} )
        endif()
    endif()
    if( atlas_TRANS_FOUND)
        if( DEFINED ENV{SABER_TEST_SPECTRALB} )
            set( SABER_TEST_SPECTRALB $ENV{SABER_TEST_SPECTRALB} )
        endif()
    endif()
    
    # Summary
    message( STATUS "SABER variables:" )
    message( STATUS "  - SABER_TEST_TIER:          ${SABER_TEST_TIER}" )
    message( STATUS "  - SABER_TEST_MPI:           ${SABER_TEST_MPI}" )
    message( STATUS "  - SABER_TEST_OMP:           ${SABER_TEST_OMP}" )
    message( STATUS "  - SABER_TEST_MULTI:         ${SABER_TEST_MULTI}" )
    message( STATUS "  - SABER_TEST_GSI:           ${SABER_TEST_GSI}" )
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
    if( SABER_TEST_MULTI )
        message( STATUS "  - TIER 1 multicores" )
        file( STRINGS testlist/saber_test_tier1-multi.txt saber_test_multi_tmp )
        list( APPEND saber_test_multi ${saber_test_multi_tmp} )
        list( APPEND saber_test_full ${saber_test_multi_tmp} )
    endif()
    if( SABER_TEST_GSI )
        message( STATUS "  - TIER 1 GSI-specific" )
        file( STRINGS testlist/saber_test_tier1-gsi.txt saber_test_gsi_tmp )
        list( APPEND saber_test_gsi ${saber_test_gsi_tmp} )
        list( APPEND saber_test_full ${saber_test_gsi_tmp} )
    endif()
    if( SABER_TEST_SPECTRALB )
        message( STATUS "  - TIER 1 SPECTRALB-specific" )
        file( STRINGS testlist/saber_test_tier1-spectralb.txt saber_test_spectralb_tmp )
        list( APPEND saber_test_spectralb ${saber_test_spectralb_tmp} )
        list( APPEND saber_test_full ${saber_test_spectralb_tmp} )
    endif()
    
    # TIER > 1
    if( SABER_TEST_TIER GREATER 1 )
        message( STATUS "  - TIER 2 base" )
        file( STRINGS testlist/saber_test_tier2.txt saber_test_tmp )
        list( APPEND saber_test ${saber_test_tmp} )
        list( APPEND saber_test_full ${saber_test_tmp} )
    endif()
    
    # Set download URL
    
    # Define branch name
    if( DEFINED ENV{SABER_TESTFILES_BRANCH} )
        # Branch specified by environment
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
            set( checksum "1")
        else()
            set( ECBUILD_DOWNLOAD_BASE_URL https://dashrepo.ucar.edu/api/v1/dataset/147b_jcsda/version/1.0.0/file )
            set( checksum "0")
        endif()
        message( STATUS "Files downloaded from: " ${ECBUILD_DOWNLOAD_BASE_URL} )
    
        # Set SABER_TESTFILES_PATH
        set( SABER_TESTFILES_PATH ${CMAKE_SOURCE_DIR}/test-data-release/saber/${GIT_TAG_FUNC} )
        message( STATUS "Files downloaded into: " ${SABER_TESTFILES_PATH} )
    
        list( APPEND SABER_DATA_DOWNLOADER_ARGS
              ${ECBUILD_DOWNLOAD_BASE_URL}
              ${CMAKE_SOURCE_DIR}/test-data-release
              "saber_testinput_tier_1.tar.gz"
              ${checksum} )
    
        # Create test-data-release in source directory
        file( MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/test-data-release )
    
        # Create download script
        set ( FILENAME saber_data_downloader.py )
        set ( SOURCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME} )
        set ( DEST_FILE ${CMAKE_BINARY_DIR}/bin/${FILENAME} )
    
        # Set configure file
        if( EXISTS "${SOURCE_FILE}.in" )
          configure_file( ${SOURCE_FILE}.in ${DEST_FILE} @ONLY )
        else()
          configure_file( ${SOURCE_FILE}    ${DEST_FILE} @ONLY )
        endif()
    
        # Add download script
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
    execute_process( COMMAND bash ${CMAKE_BINARY_DIR}/bin/saber_setup.sh
                             ${SABER_TESTFILES_PATH}/testdata
                             ${CMAKE_CURRENT_SOURCE_DIR}/testref
                             ${CMAKE_CURRENT_BINARY_DIR}
                             ${CMAKE_BINARY_DIR} )
    
    # Executables
    ecbuild_add_executable( TARGET  saber_quench_randomization.x
                            SOURCES mains/quenchRandomization.cc
                            LIBS    quench
                                    saber )
    
    ecbuild_add_executable( TARGET  saber_quench_error_covariance_training.x
                            SOURCES mains/quenchErrorCovarianceTraining.cc
                            LIBS    quench
                                    saber )
    
    ecbuild_add_executable( TARGET  saber_quench_saber_block_test.x
                            SOURCES mains/quenchSaberBlockTest.cc
                            LIBS    quench
                                    saber )
    
    ecbuild_add_executable( TARGET  saber_quench_dirac.x
                            SOURCES mains/quenchDirac.cc
                            LIBS    quench
                                    saber )
    
    # BUMP tests
    set( mpi_list 1 1 2 )
    set( omp_list 1 2 1 )
    foreach( mpixomp RANGE 2 )
        # Get mpi and omp
        list( GET mpi_list ${mpixomp} mpi )
        list( GET omp_list ${mpixomp} omp )

        # Link to yaml file
        foreach( test ${saber_test} )
            execute_process( COMMAND     sed "-e s/_MPI_/${mpi}/g;s/_OMP_/${omp}/g"
                             INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                             OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/${mpi}-${omp}/${test}.yaml )
        endforeach()

        # Randomization tests
        foreach( test ${saber_test} )
            string( FIND ${test} "quench_randomization" quench )
            if( quench MATCHES 0 )
                ecbuild_add_test( TARGET saber_test_${test}_${mpi}-${omp}
                                  MPI ${mpi}
                                  OMP ${omp}
                                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_randomization.x
                                  ARGS testinput/${mpi}-${omp}/${test}.yaml
                                  DEPENDS saber_quench_randomization.x )
            endif()
        endforeach()
        
        # Error covariance training tests
        foreach( test ${saber_test} )
            string( FIND ${test} "quench_error_covariance_training" quench )
            if( quench MATCHES 0 )
                ecbuild_add_test( TARGET saber_test_${test}_${mpi}-${omp}
                                  MPI ${mpi}
                                  OMP ${omp}
                                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_training.x
                                  ARGS testinput/${mpi}-${omp}/${test}.yaml
                                  DEPENDS saber_quench_error_covariance_training.x )
            endif()
        endforeach()
        
        # Adjoint/inverse tests
        foreach( test ${saber_test} )
            string( FIND ${test} "quench_saber_block_test" quench )
            if( quench MATCHES 0 )
                ecbuild_add_test( TARGET saber_test_${test}_${mpi}-${omp}
                                  MPI ${mpi}
                                  OMP ${omp}
                                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_saber_block_test.x
                                  ARGS testinput/${mpi}-${omp}/${test}.yaml
                                  DEPENDS saber_quench_saber_block_test.x )
            endif()
        endforeach()
        
        # Dirac tests
        foreach( test ${saber_test} )
            string( FIND ${test} "quench_dirac" quench )
            if( quench MATCHES 0 )
                ecbuild_add_test( TARGET saber_test_${test}_${mpi}-${omp}
                                  MPI ${mpi}
                                  OMP ${omp}
                                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_dirac.x
                                  ARGS testinput/${mpi}-${omp}/${test}.yaml
                                  DEPENDS saber_quench_dirac.x )
            endif()
        endforeach()
    endforeach()

    # Tests depending on GSI B Climatology
    if( SABER_TEST_GSI )
      # Tasks/threads
      set( mpi 1 )
      set( omp 2 )

      # Link to yaml files
      foreach( test ${saber_test_gsi} )
         execute_process( COMMAND     sed "-e s/_MPI_/${mpi}/g;s/_OMP_/${omp}/g"
                          INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                          OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/${mpi}-${omp}/${test}.yaml )
      endforeach()
    
      # Link to nml files
      foreach( test ${saber_test_gsi} )
          execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                           ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.nml
                           ${CMAKE_CURRENT_BINARY_DIR}/testinput/${mpi}-${omp}/${test}.nml )
      endforeach()

      # Dirac tests
      foreach( test ${saber_test_gsi} )
          string( FIND ${test} "quench_dirac" quench )
          if( quench MATCHES 0 )
              ecbuild_add_test( TARGET saber_test_${test}_${mpi}-${omp}
                                MPI ${mpi}
                                OMP ${omp}
                                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_dirac.x
                                ARGS testinput/${mpi}-${omp}/${test}.yaml
                                DEPENDS saber_quench_dirac.x )
          endif()
      endforeach()
    endif()
    
    # Tests depending on spectral transforms
    if( SABER_TEST_SPECTRALB )
      # Tasks/threads
      set( mpi 1 )
      set( omp 1 )

      # Link to yaml files
      foreach( test ${saber_test_spectralb} )
         execute_process( COMMAND     sed "-e s/_MPI_/${mpi}/g;s/_OMP_/${omp}/g"
                          INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                          OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/${mpi}-${omp}/${test}.yaml )
      endforeach()
   
      # Adjoint/inverse tests
      foreach( test ${saber_test_spectralb} )
          string( FIND ${test} "quench_saber_block_test" quench )
          if( quench MATCHES 0 )
              ecbuild_add_test( TARGET saber_test_${test}_${mpi}-${omp}
                                MPI ${mpi}
                                OMP ${omp}
                                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_saber_block_test.x
                                ARGS testinput/${mpi}-${omp}/${test}.yaml
                                DEPENDS saber_quench_saber_block_test.x )
          endif()
      endforeach()
    
      # Dirac tests
      foreach( test ${saber_test_spectralb} )
          string( FIND ${test} "quench_dirac" quench )
          if( quench MATCHES 0 )
              ecbuild_add_test( TARGET saber_test_${test}_${mpi}-${omp}
                                MPI ${mpi}
                                OMP ${omp}
                                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_dirac.x
                                ARGS testinput/${mpi}-${omp}/${test}.yaml
                                DEPENDS saber_quench_dirac.x )
          endif()
      endforeach()
    endif()
else()
    message( STATUS "OOPS missing: no testing of SABER" )
endif()
