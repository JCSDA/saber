# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

if( oops_FOUND )
    # Tests activated
    message( STATUS "SABER tests activated:" )

    # Build quench
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
    elseif( SABER_TEST_TIER EQUAL 2 )
        set( SABER_TEST_OMP 1 )
    else()
        message( FATAL_ERROR "SABER_TEST_TIER should be 1 or 2, not ${SABER_TEST_TIER}" )
    endif()

    # Default test selection variables (TIER-independent)
    set( SABER_TEST_VALGRIND 0 )
    set( SABER_TEST_GSI 0 )
    if( gsibclim_FOUND )
        set( SABER_TEST_GSI 1 )
    endif()
    set( SABER_TEST_SPECTRALB 0 )
    if( atlas_TRANS_FOUND )
        set( SABER_TEST_SPECTRALB 1 )
    endif()

    # Override test selection variables using environment variables
    if( DEFINED ENV{SABER_TEST_MPI} )
        set( SABER_TEST_MPI $ENV{SABER_TEST_MPI} )
    endif()
    if( DEFINED ENV{SABER_TEST_OMP} )
        set( SABER_TEST_OMP $ENV{SABER_TEST_OMP} )
    endif()
    if( DEFINED ENV{SABER_TEST_VALGRIND} )
        set( SABER_TEST_VALGRIND $ENV{SABER_TEST_VALGRIND} )
    endif()
    if( gsibclim_FOUND )
        if( DEFINED ENV{SABER_TEST_GSI} )
            set( SABER_TEST_GSI $ENV{SABER_TEST_GSI} )
        endif()
    endif()
    if( atlas_TRANS_FOUND )
        if( DEFINED ENV{SABER_TEST_SPECTRALB} )
            set( SABER_TEST_SPECTRALB $ENV{SABER_TEST_SPECTRALB} )
        endif()
    endif()
   
    # TIER 1
    message( STATUS "  - TIER 1 base" )
    file( STRINGS testlist/saber_test_tier1.txt saber_test )
    list( APPEND saber_test_full ${saber_test} )
    if( SABER_TEST_GSI )
        message( STATUS "  - TIER 1 GSI-specific" )
        file( STRINGS testlist/saber_test_tier1-gsi.txt saber_test )
        list( APPEND saber_test_full ${saber_test} )
    endif()
    if( SABER_TEST_SPECTRALB )
        message( STATUS "  - TIER 1 SPECTRALB-specific" )
        file( STRINGS testlist/saber_test_tier1-spectralb.txt saber_test )
        list( APPEND saber_test_full ${saber_test} )
    endif()

    # TIER > 1
    if( SABER_TEST_TIER GREATER 1 )
        message( STATUS "  - TIER 2 base" )
        file( STRINGS testlist/saber_test_tier2.txt saber_test )
        list( APPEND saber_test_full ${saber_test} )
    endif()

    # Input data
    file( STRINGS testlist/saber_data.txt saber_data )

    # Setup SABER directories and links
    message( STATUS "Setup SABER directories and links" )
    file( REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/testdata )
    file( REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/testinput )
    file( REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/testref )
    file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testdata )
    file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput )
    file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref )
    foreach( test ${saber_test_full} )
        # Make test-specific data directory
        file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testdata/${test} )

        # Link to yaml file
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                         ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${test}.yaml
                         ${CMAKE_CURRENT_BINARY_DIR}/testinput/${test}.yaml )

        # Link to reference file
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                         ${CMAKE_CURRENT_SOURCE_DIR}/testref/${test}.ref
                         ${CMAKE_CURRENT_BINARY_DIR}/testref/${test}.ref )
    endforeach()
    foreach( data ${saber_data} )
        # Link to input data file
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                         ${CMAKE_CURRENT_SOURCE_DIR}/testdata/${data}
                         ${CMAKE_CURRENT_BINARY_DIR}/testdata/${data} )
    endforeach()

    # VADER libraries
    if( vader_FOUND )
        set( vader_LIBRARIES "vader" )
    else()
        set( vader_LIBRARIES "" )
    endif()

    # Executables
    ecbuild_add_executable( TARGET  saber_quench_randomization.x
                            SOURCES mains/quenchRandomization.cc
                            LIBS    quench
                                    ${vader_LIBRARIES}
                                    saber )

    ecbuild_add_executable( TARGET  saber_quench_error_covariance_training.x
                            SOURCES mains/quenchErrorCovarianceTraining.cc
                            LIBS    quench
                                    ${vader_LIBRARIES}
                                    saber )

    ecbuild_add_executable( TARGET  saber_quench_saber_block_test.x
                            SOURCES mains/quenchSaberBlockTest.cc
                            LIBS    quench
                                    ${vader_LIBRARIES}
                                    saber )

    ecbuild_add_executable( TARGET  saber_quench_dirac.x
                            SOURCES mains/quenchDirac.cc
                            LIBS    quench
                                    ${vader_LIBRARIES}
                                    saber )

    ecbuild_add_executable( TARGET  saber_quench_convertstate.x
                            SOURCES mains/quenchConvertState.cc
                            LIBS    quench
                                    ${vader_LIBRARIES}
                                    saber )

    # List of MPI/OpenMP configurations to test
    list( APPEND mpi_list 1)
    list( APPEND omp_list 1)
    if( SABER_TEST_MPI)
        list( APPEND mpi_list 2)
        list( APPEND omp_list 1)
    endif()
    if( SABER_TEST_OMP)
        list( APPEND mpi_list 1)
        list( APPEND omp_list 2)
    endif()

    # Loop over MPI/OpenMP configurations
    foreach( mpi omp IN ZIP_LISTS mpi_list omp_list )
        message( STATUS "Tests parallel configuration: MPI: ${mpi} / OMP: ${omp}" )
 
        # Randomization tests
        foreach( test ${saber_test_full} )
            string( FIND ${test} "randomization" result )
            if( result MATCHES 0 )
                ecbuild_add_test( TARGET saber_test_${test}_${mpi}-${omp}
                                  MPI ${mpi}
                                  OMP ${omp}
                                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_randomization.x
                                  ARGS testinput/${test}.yaml
                                  DEPENDS saber_quench_randomization.x )
            endif()
        endforeach()
    
        # Error covariance training tests
        foreach( test ${saber_test_full} )
            string( FIND ${test} "error_covariance_training" result )
            if( result MATCHES 0 )
                ecbuild_add_test( TARGET saber_test_${test}_${mpi}-${omp}
                                  MPI ${mpi}
                                  OMP ${omp}
                                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_training.x
                                  ARGS testinput/${test}.yaml
                                  DEPENDS saber_quench_error_covariance_training.x )
            endif()
        endforeach()
    
        # Adjoint/inverse tests
        foreach( test ${saber_test_full} )
            string( FIND ${test} "saber_block_test" result )
            if( result MATCHES 0 )
                ecbuild_add_test( TARGET saber_test_${test}_${mpi}-${omp}
                                  MPI ${mpi}
                                  OMP ${omp}
                                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_saber_block_test.x
                                  ARGS testinput/${test}.yaml
                                  DEPENDS saber_quench_saber_block_test.x )
            endif()
        endforeach()
    
        # Dirac tests
        foreach( test ${saber_test_full} )
            string( FIND ${test} "dirac" result )
            if( result MATCHES 0 )
                ecbuild_add_test( TARGET saber_test_${test}_${mpi}-${omp}
                                  MPI ${mpi}
                                  OMP ${omp}
                                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_dirac.x
                                  ARGS testinput/${test}.yaml
                                  DEPENDS saber_quench_dirac.x )
            endif()
        endforeach()

        # ConvertState tests
        foreach( test ${saber_test_full} )
            string( FIND ${test} "convertstate" result )
            if( result MATCHES 0 )
                ecbuild_add_test( TARGET saber_test_${test}_${mpi}-${omp}
                                  MPI ${mpi}
                                  OMP ${omp}
                                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_dirac.x
                                  ARGS testinput/${test}.yaml
                                  DEPENDS saber_quench_dirac.x )
            endif()
        endforeach()

    endforeach()
else()
    # Tests not activated
    message( STATUS "SABER tests not activated" )
endif()
