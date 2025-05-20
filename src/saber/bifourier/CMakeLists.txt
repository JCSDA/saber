# (C) Copyright 2025 Meteorologisk Institutt.
#

list(APPEND bifourier_src_files_list
# Biperiodization block
Biperiodization.h
Biperiodization.cc
BiperiodizationImpl.h
BiperiodizationImpl.cc

# Reduced wind to geographical wind
RedWindToGeoWind.cc
RedWindToGeoWind.h
)

if( FFTW_FOUND )
    list(APPEND bifourier_src_files_list
    # Bi-Fourier AROME legacy file format reading
    bifourier_arome_legacy.h
    BifourierAromeLegacy.cc
    BifourierAromeLegacy.h
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

    # Bi-Fourier temperature/surface pressure split block
    BifourierSplitTPs.cc
    BifourierSplitTPs.h

    # Bi-Fourier spectral transform object
    BifourierTransform.cc
    BifourierTransform.h

    # Bi-Fourier spectral transforms store
    BifourierTransformStore.cc
    BifourierTransformStore.h

    # Bi-Fourier utilities
    BifourierUtilities.cc
    BifourierUtilities.h

    # Bi-Fourier vorticity/divergence to reduced wind
    BifourierVorDivToRedWind.cc
    BifourierVorDivToRedWind.h

    # Bi-Fourier vorticity to balance pressure
    BifourierVorToPb.cc
    BifourierVorToPb.h
    )
endif()

set( bifourier_src_files ${bifourier_src_files_list} PARENT_SCOPE )
