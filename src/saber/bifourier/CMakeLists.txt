# (C) Copyright 2025 Meteorologisk Institutt.
#

list(APPEND bifourier_src_files_list
# Biperiodization block
Biperiodization.h
Biperiodization.cc
BiperiodizationImpl.h
BiperiodizationImpl.cc
)

if( FFTW_FOUND OR ( ECTRANS_FOUND AND ectrans_HAVE_ETRANS AND ectrans_HAVE_TRANSI ) )
    list(APPEND bifourier_src_files_list
    # Bi-Fourier AROME balance block
    BifourierAromeBalance.cc
    BifourierAromeBalance.h

    # Bi-Fourier AROME covariance block
    BifourierAromeCovariance.cc
    BifourierAromeCovariance.h

    # Bi-Fourier AROME legacy file format reading
    bifourier_arome_legacy.h
    bifourier_arome_legacy_interface.F90
    bifourier_arome_legacy_mod.F90

    # Bi-Fourier balance block
    BifourierBalance.cc
    BifourierBalance.h

    # Bi-Fourier covariance block
    BifourierCovariance.cc
    BifourierCovariance.h

    # Bi-Fourier grid to spectral block
    BifourierGridToSpectral.cc
    BifourierGridToSpectral.h

    # Bi-Fourier identity block (in spectral space)
    BifourierID.cc
    BifourierID.h

    # Bi-Fourier spectral to grid block
    BifourierSpectralToGrid.cc
    BifourierSpectralToGrid.h

    # Bi-Fourier spectral vorticity/divergence to grid wind
    BifourierSpectralVorDivToGridWind.cc
    BifourierSpectralVorDivToGridWind.h

    # Bi-Fourier spectral transform object, base class
    BifourierTransformBase.cc
    BifourierTransformBase.h

    # Bi-Fourier spectral transforms store
    BifourierTransformStore.cc
    BifourierTransformStore.h

    # Bi-Fourier utilities
    BifourierUtilities.cc
    BifourierUtilities.h
    )
endif()

# Add FFT backends
message(STATUS "Available FFT backends in BifourierTransform:")
if( FFTW_FOUND )
    # FFTW backend
    message(STATUS "- fftw")
    list(APPEND bifourier_src_files_list
    BifourierTransformFFTW.cc
    BifourierTransformFFTW.h
    )
endif()
if( ECTRANS_FOUND AND ectrans_HAVE_ETRANS AND ectrans_HAVE_TRANSI )
    # ECTRANS backend
    message(STATUS "- ectrans")
    list(APPEND bifourier_src_files_list
    BifourierTransformECTRANS.cc
    BifourierTransformECTRANS.h
    )
endif()

set( bifourier_src_files ${bifourier_src_files_list} PARENT_SCOPE )
