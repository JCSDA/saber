# (C) Copyright 2017-2023 UCAR.
# (C) Copyright 2023 Meteorologisk Institutt
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Tests activated
message( STATUS "SABER tests activated:" )

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
set( SABER_TEST_FASTLAM 1 )
set( SABER_TEST_GSI_GEOS 0 )
if( gsibec_FOUND )
    set( SABER_TEST_GSI_GEOS 1 )
endif()
set( SABER_TEST_GSI_GFS 0 )
if( gsibec_FOUND AND sp_FOUND )
    set( SABER_TEST_GSI_GFS 1 )
endif()
set( SABER_TEST_SPECTRALB 0 )
if( atlas_TRANS_FOUND OR atlas_ECTRANS_FOUND )
    set( SABER_TEST_SPECTRALB 1 )
endif()
set( SABER_TEST_VADER 0 )
if( vader_FOUND )
    set( SABER_TEST_VADER 1 )
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
if( DEFINED ENV{SABER_TEST_FASTLAM} )
    set( SABER_TEST_FASTLAM $ENV{SABER_TEST_FASTLAM} )
endif()
if( gsibec_FOUND )
    if( DEFINED ENV{SABER_TEST_GSI_GEOS} )
        set( SABER_TEST_GSI_GEOS $ENV{SABER_TEST_GSI_GEOS} )
    endif()
    if( sp_FOUND )
        if( DEFINED ENV{SABER_TEST_GSI_GFS} )
            set( SABER_TEST_GSI_GFS $ENV{SABER_TEST_GSI_GFS} )
        endif()
    endif()
endif()
if( atlas_TRANS_FOUND OR atlas_ECTRANS_FOUND )
    if( DEFINED ENV{SABER_TEST_SPECTRALB} )
        set( SABER_TEST_SPECTRALB $ENV{SABER_TEST_SPECTRALB} )
    endif()
endif()
if( vader_FOUND )
    if( DEFINED ENV{SABER_TEST_VADER} )
        set( SABER_TEST_VADER $ENV{SABER_TEST_VADER} )
    endif()
endif()

# TIER 1
message( STATUS "  - TIER 1 base" )
file( STRINGS testlist/saber_test_tier1.txt saber_test )
list( APPEND saber_test_full ${saber_test} )
if( SABER_TEST_FASTLAM )
    message( STATUS "  - TIER 1 FastLAM-specific" )
    file( STRINGS testlist/saber_test_tier1-fastlam.txt saber_test )
    list( APPEND saber_test_full ${saber_test} )
    if( FFTW_FOUND )
        message( STATUS "  - TIER 1 FastLAM-specific - FFTW" )
        file( STRINGS testlist/saber_test_tier1-fastlam-fftw.txt saber_test )
        list( APPEND saber_test_full ${saber_test} )
    endif()
endif()
if( SABER_TEST_GSI_GEOS )
    message( STATUS "  - TIER 1 GSI-GEOS-specific" )
    file( STRINGS testlist/saber_test_tier1-gsi-geos.txt saber_test )
    list( APPEND saber_test_full ${saber_test} )
endif()
if( SABER_TEST_GSI_GFS )
    message( STATUS "  - TIER 1 GSI-GFS-specific" )
    file( STRINGS testlist/saber_test_tier1-gsi-gfs.txt saber_test )
    list( APPEND saber_test_full ${saber_test} )
endif()
if( SABER_TEST_SPECTRALB )
    message( STATUS "  - TIER 1 SPECTRALB-specific" )
    file( STRINGS testlist/saber_test_tier1-spectralb.txt saber_test )
    list( APPEND saber_test_full ${saber_test} )
endif()
if( SABER_TEST_VADER )
    message( STATUS "  - TIER 1 VADER-specific" )
    file( STRINGS testlist/saber_test_tier1_vader.txt saber_test )
    list( APPEND saber_test_full ${saber_test} )
    if( atlas_TRANS_FOUND OR atlas_ECTRANS_FOUND )
       file( STRINGS testlist/saber_test_tier1_vader_with_trans.txt saber_test )
       list( APPEND saber_test_full ${saber_test} )
    endif()
    if ( SABER_TEST_SPECTRALB )
       file( STRINGS testlist/saber_test_tier1-spectralb-vader.txt saber_test )
       list( APPEND saber_test_full ${saber_test} )
    endif()
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
    if( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/testref/${test}.ref )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                         ${CMAKE_CURRENT_SOURCE_DIR}/testref/${test}.ref
                         ${CMAKE_CURRENT_BINARY_DIR}/testref/${test}.ref )
    endif()
endforeach()
foreach( data ${saber_data} )
    # Link to input data file
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                     ${CMAKE_CURRENT_SOURCE_DIR}/testdata/${data}
                     ${CMAKE_CURRENT_BINARY_DIR}/testdata/${data} )
endforeach()

# List of MPI/OpenMP configurations to test
list( APPEND mpi_list 1)
list( APPEND omp_list 1)
if( SABER_TEST_MPI )
    list( APPEND mpi_list 2)
    list( APPEND mpi_list 4)
    list( APPEND omp_list 1)
    list( APPEND omp_list 1)
endif()
if( SABER_TEST_OMP )
    list( APPEND mpi_list 1)
    list( APPEND omp_list 2)
endif()

# Executables list
list( APPEND exe_list convertstate randomization error_covariance_training process_perts dirac )

# Loop over MPI/OpenMP configurations
foreach( mpi omp IN ZIP_LISTS mpi_list omp_list )
    message( STATUS "Tests parallel configuration: MPI: ${mpi} / OMP: ${omp}" )

    # Loop over executables
    foreach( exe ${exe_list} )
        if( "${exe}" STREQUAL "convertstate" )
            set( exename "convertstate" )
        elseif ("${exe}" STREQUAL "process_perts" )
            set( exename "process_perts" )
        else()
            set( exename "error_covariance_toolbox" )
        endif()
        foreach( test ${saber_test_full} )
            string( FIND ${test} ${exe} result )
            if( result MATCHES 0 )
                # Special check for doc and tutorial, only run with 1 MPI / 1 OMP
                string( FIND ${test} "doc_" doc_result )
                string( FIND ${test} "tutorial_" tutorial_result )
                if( ( doc_result EQUAL -1 AND tutorial_result EQUAL -1 ) OR ( ${mpi} EQUAL 1 AND ${omp} EQUAL 1 ) )
                    set( docTutorialCheck true CACHE BOOL "" FORCE )
                else()
                    set( docTutorialCheck false CACHE BOOL "" FORCE )
                endif()

                # Special check for GSI GFS, only run with 1 MPI / 1 OMP (TODO: should be removed at some point)
                string( FIND ${test} "gsi_gfs" gsi_gfs_result )
                if( gsi_gfs_result EQUAL -1 OR ( ${mpi} EQUAL 1 AND ${omp} EQUAL 1 ) )
                    set( gsiGfsCheck true CACHE BOOL "" FORCE )
                else()
                    set( gsiGfsCheck false CACHE BOOL "" FORCE )
                endif()

                # Special check for 4D tests
                string( FIND ${test} "_4d" fourd_result )
                if( fourd_result EQUAL -1 OR ( fourd_result GREATER_EQUAL 0 AND ${mpi} GREATER 1 ) )
                    set( 4dCheck true CACHE BOOL "" FORCE )
                else()
                    set( 4dCheck false CACHE BOOL "" FORCE )
                endif()

                # Special check for parallel_hybrid tests, only run with 2 MPI and more
                # Special check for non-parallel hybrid tests, only run with 2 MPI and less
                string( FIND ${test} "parallel_hybrid" par_hyb_result )
                if ( par_hyb_result EQUAL -1 AND ${mpi} LESS_EQUAL 2 )
                  set (parallelHybridCheck true CACHE BOOL "" FORCE )
                else()
                  if ( par_hyb_result GREATER_EQUAL 0 AND ${mpi} GREATER_EQUAL 2 )
                    set (parallelHybridCheck true CACHE BOOL "" FORCE )
                  else()
                    set (parallelHybridCheck false CACHE BOOL "" FORCE )
                  endif()
                endif()

                if( docTutorialCheck AND gsiGfsCheck AND 4dCheck AND parallelHybridCheck )
                    # Get dependencies
                    file( STRINGS testdeps/${test}.txt deps )
                    set( deps_list "" )
                    list( APPEND deps_list ${deps} )
                    list( LENGTH deps_list deps_length )
                    if( ${deps_length} GREATER 0 )
                        list( TRANSFORM deps_list PREPEND saber_test_ )
                        list( TRANSFORM deps_list APPEND _${mpi}-${omp} )
                        if( ${mpi} GREATER 1 )
                            set( deps_list_1 "" )
                            list( APPEND deps_list_1 ${deps} )
                            list( TRANSFORM deps_list_1 PREPEND saber_test_ )
                            list( TRANSFORM deps_list_1 APPEND _1-1 )
                            list( APPEND deps_list ${deps_list_1} )
                        endif()
                    endif()

                    # Add test
                    ecbuild_add_test( TARGET saber_test_${test}_${mpi}-${omp}
                                      MPI ${mpi}
                                      OMP ${omp}
                                      COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_${exename}.x
                                      ARGS testinput/${test}.yaml
                                      DEPENDS saber_quench_${exename}.x
                                      TEST_DEPENDS ${deps_list} )
                endif()
            endif()
        endforeach()
    endforeach()

    # GSI scaling optimization for HDIAG
    foreach( test ${saber_test_full} )
        string( FIND ${test} "optimization_bump_hdiag_gsi" result )
        if( result MATCHES 0 AND ${mpi} LESS_EQUAL 2 )
            # Get dependencies
            file( STRINGS testdeps/${test}.txt deps )
            set( deps_list "" )
            list( APPEND deps_list ${deps} )
            list( LENGTH deps_list deps_length )
            if( ${deps_length} GREATER 0 )
                list( TRANSFORM deps_list PREPEND saber_test_ )
                list( TRANSFORM deps_list APPEND _${mpi}-${omp} )
            endif()

            # Add test
            ecbuild_add_test( TARGET saber_test_${test}_${mpi}-${omp}
                              MPI ${mpi}
                              OMP ${omp}
                              COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                              ARGS testinput/${yaml_prefix}${test}.yaml
                              DEPENDS saber_quench_error_covariance_toolbox.x
                              TEST_DEPENDS ${deps_list} )
        endif()
    endforeach()

    # Compare diagnostics
    foreach( test ${saber_test_full} )
        string( FIND ${test} "compare_diagnostics" result )
        if( result MATCHES 0 )
            if ( ${mpi} EQUAL 1 AND ${omp} EQUAL 1 )
                # Add test
                ecbuild_add_test( TARGET saber_test_${test}
                                  TYPE SCRIPT
                                  COMMAND ${CMAKE_BINARY_DIR}/bin/saber_compare_dirac_diagnostics.py
                                  ARGS testinput/${test}.yaml )
            endif()
        endif()
    endforeach()
endforeach()

# BUMP interpolator test
set( mpi 2 )
set( omp 1 )
add_fctest( TARGET          saber_test_interpolatorbump
            MPI             ${mpi}
            OMP             ${omp}
            LINKER_LANGUAGE Fortran
            SOURCES         fctest/fctest_interpolatorbump.F90
            LIBS            saber )

# NICAS square-root test
set( mpi 2 )
set( omp 1 )
add_fctest( TARGET          saber_test_nicas_sqrt
            MPI             ${mpi}
            OMP             ${omp}
            LINKER_LANGUAGE Fortran
            SOURCES         fctest/fctest_nicas_sqrt.F90
            LIBS            saber )
