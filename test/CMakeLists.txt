# (C) Copyright 2017-2023 UCAR.
# (C) Copyright 2023 Meteorologisk Institutt
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

if( ENABLE_QUENCH )
  # Tests activated
  message( STATUS "SABER-QUENCH tests activated:" )

  # Use find_branch function in jedi-cmake
  include( git_functions )

  # Default test selection variables (TIER-independent)
  set( SABER_TEST_FASTLAM 1 )
  set( SABER_TEST_BIFOURIER 0 )
  if ( FFTW_FOUND OR (ECTRANS_FOUND AND ectrans_HAVE_ETRANS AND ectrans_HAVE_TRANSI ) )
    set( SABER_TEST_BIFOURIER 1 )
  endif()
  set( SABER_TEST_BUMP 0 )
  if( ENABLE_BUMP )
      set( SABER_TEST_BUMP 1 )
  endif()
  set( SABER_TEST_GSI_GEOS 0 )
  if( gsibec_FOUND )
      set( SABER_TEST_GSI_GEOS 1 )
  endif()
  set( SABER_TEST_GSI_GFS 0 )
  if( gsibec_FOUND AND ip_FOUND )
      set( SABER_TEST_GSI_GFS 1 )
  endif()
  set( SABER_TEST_SPECTRALB 0 )
  if( atlas_TRANS_FOUND OR atlas_ECTRANS_FOUND )
      set( SABER_TEST_SPECTRALB 1 )
  endif()
  set( SABER_TEST_VADER 1 )
  set( SABER_TEST_TORCHBALANCE 0 )
  if( ENABLE_TORCH )
      set( SABER_TEST_TORCHBALANCE 1 )
  endif()

  # Override test selection variables using environment variables
  if( DEFINED ENV{SABER_TEST_FASTLAM} )
      set( SABER_TEST_FASTLAM $ENV{SABER_TEST_FASTLAM} )
  endif()
  if( DEFINED ENV{SABER_TEST_BIFOURIER} )
      set( SABER_TEST_BIFOURIER $ENV{SABER_TEST_BIFOURIER} )
  endif()
  if( ENABLE_BUMP )
      if( DEFINED ENV{SABER_TEST_BUMP} )
          set( SABER_TEST_BUMP $ENV{SABER_TEST_BUMP} )
      endif()
  endif()
  if( gsibec_FOUND )
      if( DEFINED ENV{SABER_TEST_GSI_GEOS} )
          set( SABER_TEST_GSI_GEOS $ENV{SABER_TEST_GSI_GEOS} )
      endif()
      if( ip_FOUND )
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
  if( DEFINED ENV{SABER_TEST_VADER} )
      set( SABER_TEST_VADER $ENV{SABER_TEST_VADER} )
  endif()

  # BUMP needed for FASTLAM tests
  if( NOT SABER_TEST_BUMP )
      set( SABER_TEST_FASTLAM 0 )
  endif()

  if( SABER_TEST_TORCHBALANCE )
      # Add test to generate emulator (runs before torchbalance tests)
      ecbuild_add_test( TARGET generate_torchbalance_test_emulator
                        TYPE SCRIPT
                        COMMAND ${VENV_PYTHON}
                        ARGS ${CMAKE_SOURCE_DIR}/saber/src/saber/torchbalance/create_test_emulator.py
                             ${CMAKE_CURRENT_BINARY_DIR}/testdata/dirac_torchbalance )
  endif()

  # Setup SABER top-level directories
  file( REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/testdata )
  file( REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/testinput )
  file( REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/testref )
  file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testdata )
  file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput )
  file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref )

  # Wrapper around ecbuild_add_test that also creates the per-test testdata/
  # directory and symlinks the yaml input and reference files into testinput/
  # and testref/. The test name is derived from TARGET by stripping the
  # "saber_" prefix and the trailing "_<mpi>-<omp>" suffix.
  macro( saber_add_test )
    cmake_parse_arguments( _SAT "" "TARGET" "" ${ARGN} )

    # Extract test name: strip saber_ prefix and optional _N-N config suffix
    string( REGEX REPLACE "^saber_" "" _sat_name "${_SAT_TARGET}" )
    string( REGEX REPLACE "_[0-9]+-[0-9]+$" "" _sat_name "${_sat_name}" )

    # Only perform filesystem setup once per unique test name
    if( NOT DEFINED _SABER_SETUP_DONE_${_sat_name} )

      set( _SABER_SETUP_DONE_${_sat_name} TRUE )
      # Create test-specific data directory
      file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testdata/${_sat_name} )

      # Symlink yaml input file
      file( CREATE_LINK
            ${CMAKE_CURRENT_SOURCE_DIR}/testinput/${_sat_name}.yaml
            ${CMAKE_CURRENT_BINARY_DIR}/testinput/${_sat_name}.yaml SYMBOLIC )

      # Symlink reference file if it exists
      if( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/testref/${_sat_name}.ref )
        file( CREATE_LINK
              ${CMAKE_CURRENT_SOURCE_DIR}/testref/${_sat_name}.ref
              ${CMAKE_CURRENT_BINARY_DIR}/testref/${_sat_name}.ref SYMBOLIC )
      endif()
    endif()

    ecbuild_add_test( ${ARGN} )
  endmacro()

  # Path to the saber test files
  set (jedi_model_data_tier_1_saber "${CMAKE_SOURCE_DIR}/jedi-model-data/testinput_tier_1/saber")

  # Input data
  file( STRINGS testlist/saber_data.txt saber_data )

  foreach( data ${saber_data} )
      # Link to input data file
      file( CREATE_LINK
            ${jedi_model_data_tier_1_saber}/${data}
            ${CMAKE_CURRENT_BINARY_DIR}/testdata/${data} SYMBOLIC )
  endforeach()

  # --- Base tests ---
  include( ${CMAKE_CURRENT_LIST_DIR}/testlist/saber_base.cmake )

  # --- Coupled tests ---
  include( ${CMAKE_CURRENT_LIST_DIR}/testlist/saber_coupled.cmake)

  # --- Fast LAM tests ---
  if( SABER_TEST_FASTLAM )
    include( ${CMAKE_CURRENT_LIST_DIR}/testlist/saber_fastlam.cmake )
    include( ${CMAKE_CURRENT_LIST_DIR}/testlist/saber_fastlam-regional_interpolation.cmake )
    if( FFTW_FOUND )
      include( ${CMAKE_CURRENT_LIST_DIR}/testlist/saber_fastlam-fftw.cmake )
    endif()
  endif()

  # --- Bifourier tests ---
  if( SABER_TEST_BIFOURIER )
    include( ${CMAKE_CURRENT_LIST_DIR}/testlist/saber_bifourier.cmake )
    if( ECTRANS_FOUND AND ectrans_HAVE_ETRANS AND ectrans_HAVE_TRANSI )
      include( ${CMAKE_CURRENT_LIST_DIR}/testlist/saber_bifourier-ectrans.cmake )
    endif()
  endif()

  # --- BUMP tests ---
  if( SABER_TEST_BUMP )
    include( ${CMAKE_CURRENT_LIST_DIR}/testlist/saber_bump.cmake )
  endif()

  # --- GSI tests ---
  if( SABER_TEST_GSI_GEOS )
    include( ${CMAKE_CURRENT_LIST_DIR}/testlist/saber_gsi-geos.cmake )
  endif()
  if( SABER_TEST_GSI_GFS )
    include( ${CMAKE_CURRENT_LIST_DIR}/testlist/saber_gsi-gfs.cmake )
  endif()

  # --- SpectralB tests ---
  if( SABER_TEST_SPECTRALB )
    include( ${CMAKE_CURRENT_LIST_DIR}/testlist/saber_spectralb.cmake )
  endif()

  # --- VADER tests ---
  if( SABER_TEST_VADER )
    include( ${CMAKE_CURRENT_LIST_DIR}/testlist/saber_vader.cmake )
    if( atlas_TRANS_FOUND OR atlas_ECTRANS_FOUND )
      include( ${CMAKE_CURRENT_LIST_DIR}/testlist/saber_vader_with_trans.cmake )
    endif()
    if( SABER_TEST_SPECTRALB )
      include( ${CMAKE_CURRENT_LIST_DIR}/testlist/saber_spectralb-vader.cmake )
    endif()
  endif()

  # --- TorchBalance tests ---
  if( SABER_TEST_TORCHBALANCE )
    include( ${CMAKE_CURRENT_LIST_DIR}/testlist/saber_torchbalance.cmake )
  endif()
else()
  # Tests deactivated
  message( STATUS "SABER-QUENCH tests deactivated" )
endif()

if( SABER_TEST_BUMP )
    # NICAS square-root test
    add_fctest( TARGET          saber_nicas_sqrt
                MPI             2
                OMP             1
                LABELS          tier2
                LINKER_LANGUAGE ${FORTRAN_LINKER_LANGUAGE}
                SOURCES         fctest/fctest_nicas_sqrt.F90
                LIBS            saber )
endif()
