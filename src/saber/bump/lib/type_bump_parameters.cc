/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/lib/type_bump_parameters.h"

#include <csignal>
#include <iostream>

#include "eckit/config/LocalConfiguration.h"

namespace bump_lib {

// -----------------------------------------------------------------------------

template <class T>
void param(std::pair<const char*, T> pair, eckit::LocalConfiguration & config) {
  config.set(pair.first, pair.second);
}

// -----------------------------------------------------------------------------

void bump_config_init_f90(eckit::LocalConfiguration * config) {
  //
  if (!config->empty()) {
    std::cout << "bump_config_init_f90: LocalConfiguration should be empty" << std::endl;
    std::abort();
  }

  // General section
  GeneralDef generalDef;
  eckit::LocalConfiguration generalConf;
  // Add colors to the log (for display on terminal)
  param(generalDef.color_log, generalConf);
  // Stream test messages into a dedicated channel
  param(generalDef.testing, generalConf);
  // Default seed for random numbers
  param(generalDef.default_seed, generalConf);
  // Inter-compilers reproducibility
  param(generalDef.repro_ops, generalConf);
  // Reproducibility threshold
  param(generalDef.repro_th, generalConf);
  // Universe radius [in meters]
  param(generalDef.universe_radius, generalConf);
  // Sampling method
  param(generalDef.sampling_method, generalConf);

  // I/O section
  IODef ioDef;
  eckit::LocalConfiguration ioConf;
  // Data directory
  param(ioDef.data_directory, ioConf);
  // Data prefix
  param(ioDef.files_prefix, ioConf);
  // Write in new files
  param(ioDef.new_files, ioConf);
  // Parallel NetCDF I/O
  param(ioDef.parallel_netcdf, ioConf);
  // Number of I/O processors
  param(ioDef.nprocio, ioConf);
  // Sampling file
  param(ioDef.fname_samp, ioConf);
  // Vertical balance file
  param(ioDef.fname_vbal, ioConf);
  // NICAS file
  param(ioDef.fname_nicas, ioConf);
  // Psichitouv transform file
  param(ioDef.fname_wind, ioConf);
  // GSI data file
  param(ioDef.fname_gsi_data, ioConf);
  // GSI namelist
  param(ioDef.fname_gsi_nam, ioConf);

  // Drivers section
  DriversDef driversDef;
  eckit::LocalConfiguration driversConf;
  // Compute covariance, ensemble 1
  param(driversDef.compute_cov1, driversConf);
  // Compute covariance, ensemble 2
  param(driversDef.compute_cov2, driversConf);
  // Compute correlation, ensemble 1
  param(driversDef.compute_cor1, driversConf);
  // Compute correlation, ensemble 2
  param(driversDef.compute_cor2, driversConf);
  // Compute localization, ensemble 1
  param(driversDef.compute_loc1, driversConf);
  // Compute localization, ensemble 2
  param(driversDef.compute_loc2, driversConf);
  // Compute hybrid weights
  param(driversDef.compute_hyb, driversConf);
  // Hybrid term source ('randomized static' or 'lowres ensemble')
  param(driversDef.hybrid_source, driversConf);
  // Multivariate strategy ('univariate', 'duplicated', 'duplicated and weighted' or 'crossed')
  param(driversDef.strategy, driversConf);
  // New normality test
  param(driversDef.new_normality, driversConf);
  // Read local sampling
  param(driversDef.load_samp_local, driversConf);
  // Read global sampling
  param(driversDef.load_samp_global, driversConf);
  // Write local sampling
  param(driversDef.write_samp_local, driversConf);
  // Write global sampling
  param(driversDef.write_samp_global, driversConf);
  // Write sampling grids
  param(driversDef.write_samp_grids, driversConf);
  // New vertical covariance
  param(driversDef.new_vbal_cov, driversConf);
  // Read local vertical covariance
  param(driversDef.load_vbal_cov, driversConf);
  // Write local vertical covariancee
  param(driversDef.write_vbal_cov, driversConf);
  // Compute vertical balance operator
  param(driversDef.new_vbal, driversConf);
  // Read local vertical balance operator
  param(driversDef.load_vbal, driversConf);
  // Write vertical balance operator
  param(driversDef.write_vbal, driversConf);
  // Compute variance
  param(driversDef.new_var, driversConf);
  // Compute moments
  param(driversDef.new_mom, driversConf);
  // Read sampling moments
  param(driversDef.load_mom, driversConf);
  // Write sampling moments
  param(driversDef.write_mom, driversConf);
  // Write HDIAG diagnostics
  param(driversDef.write_hdiag, driversConf);
  // Write HDIAG components detail
  param(driversDef.write_hdiag_detail, driversConf);
  // Read universe radius
  param(driversDef.load_universe_radius, driversConf);
  // Write universe radius
  param(driversDef.write_universe_radius, driversConf);
  // Compute NICAS
  param(driversDef.new_nicas, driversConf);
  // Read local NICAS parameters
  param(driversDef.load_nicas_local, driversConf);
  // Read global NICAS parameters
  param(driversDef.load_nicas_global, driversConf);
  // Write local NICAS parameters
  param(driversDef.write_nicas_local, driversConf);
  // Write global NICAS parameters
  param(driversDef.write_nicas_global, driversConf);
  // Write NICAS grids
  param(driversDef.write_nicas_grids, driversConf);
  // Write NICAS steps
  param(driversDef.write_nicas_steps, driversConf);
  // Compute wind transform
  param(driversDef.new_wind, driversConf);
  // Read local wind transform
  param(driversDef.load_wind_local, driversConf);
  // Write local wind transform
  param(driversDef.write_wind_local, driversConf);
  // Test vertical balance inverse
  param(driversDef.check_vbal, driversConf);
  // Test adjoints
  param(driversDef.check_adjoints, driversConf);
  // Test NICAS normalization (number of tests)
  param(driversDef.check_normalization, driversConf);
  // Test NICAS application on diracs
  param(driversDef.check_dirac, driversConf);
  // Test NICAS randomization
  param(driversDef.check_randomization, driversConf);
  // Test HDIAG-NICAS consistency
  param(driversDef.check_consistency, driversConf);
  // Test HDIAG optimality
  param(driversDef.check_optimality, driversConf);
  // Interpolate vertical balance, standard-deviation or length-scales from GSI data
  param(driversDef.from_gsi, driversConf);

  // Model section
  ModelDef modelDef;
  eckit::LocalConfiguration modelConf;
  // Level for 2D variables ('first' or 'last')
  param(modelDef.lev2d, modelConf);
  // Check that sampling couples and interpolations do not cross mask boundaries
  param(modelDef.mask_check, modelConf);

  // Ensemble sizes section
  EnsembleSizesDef ensembleSizesDef;
  eckit::LocalConfiguration ensembleSizesConf;
  // Ensemble 1 size
  param(ensembleSizesDef.ens1_ne, ensembleSizesConf);
  // Ensemble 1 sub-ensembles number
  param(ensembleSizesDef.ens1_nsub, ensembleSizesConf);
  // Ensemble 2 size
  param(ensembleSizesDef.ens2_ne, ensembleSizesConf);
  // Ensemble 2 sub-ensembles number
  param(ensembleSizesDef.ens2_nsub, ensembleSizesConf);

  // Sampling section
  SamplingDef samplingDef;
  eckit::LocalConfiguration samplingConf;
  // Computation grid size
  param(samplingDef.nc1, samplingConf);
  // Diagnostic grid size
  param(samplingDef.nc2, samplingConf);
  // Number of distance classes
  param(samplingDef.nc3, samplingConf);
  // Number of angular sectors
  param(samplingDef.nc4, samplingConf);
  // Class size (for sam_type='hor'), should be larger than the typical grid cell size [in meters]
  param(samplingDef.dc, samplingConf);
  // Reduced number of levels for diagnostics
  param(samplingDef.nl0r, samplingConf);
  // Activate local diagnostics
  param(samplingDef.local_diag, samplingConf);
  // Local diagnostics calculation radius [in meters]
  param(samplingDef.local_rad, samplingConf);
  // Local diagnostics calculation latitude band half-width [in degrees]
  param(samplingDef.local_dlat, samplingConf);
  // Diagnostic draw type ('random' or 'octahedral')
  param(samplingDef.draw_type, samplingConf);
  // Maximum number of random number draws
  param(samplingDef.irmax, samplingConf);
  // Vertical balance C2B to C0A interpolation type ('c0': C0 mesh-based, 'c1': C1 mesh-based
  // or 'si': smooth interpolation)
  param(samplingDef.interp_type, samplingConf);
  // Threshold on vertically contiguous points for the mask (0 to skip the test)
  param(samplingDef.ncontig_th, samplingConf);

  // Diagnostics section
  DiagnosticsDef diagnosticsDef;
  eckit::LocalConfiguration diagnosticsConf;
  // Ensemble size
  param(diagnosticsDef.ne, diagnosticsConf);
  // Ensemble size of the hybrid term
  param(diagnosticsDef.ne_lr, diagnosticsConf);
  // Gaussian approximation for asymptotic quantities
  param(diagnosticsDef.gau_approx, diagnosticsConf);
  // Threshold on generalized kurtosis (3.0 = Gaussian distribution)
  param(diagnosticsDef.gen_kurt_th, diagnosticsConf);
  // Number of bins for averaged statistics histograms
  param(diagnosticsDef.avg_nbins, diagnosticsConf);
  // Support radius scaling in CMAT from HDIAG
  param(diagnosticsDef.lengths_scaling, diagnosticsConf);

  // Vertical balance section
  VerticalBalanceDef verticalBalanceDef;
  eckit::LocalConfiguration verticalBalanceConf;
  // Pseudo-inverse for auto-covariance
  param(verticalBalanceDef.vbal_pseudo_inv, verticalBalanceConf);
  // Dominant mode for pseudo-inverse
  param(verticalBalanceDef.vbal_pseudo_inv_mmax, verticalBalanceConf);
  // Variance threshold to compute the dominant mode for pseudo-inverse
  param(verticalBalanceDef.vbal_pseudo_inv_var_th, verticalBalanceConf);
  // Identity vertical balance for tests
  param(verticalBalanceDef.vbal_id, verticalBalanceConf);

  // Variance section
  VarianceDef varianceDef;
  eckit::LocalConfiguration varianceConf;
  // Force specific variance
  param(varianceDef.forced_var, varianceConf);
  // Filter variance
  param(varianceDef.var_filter, varianceConf);
  // Number of iterations for the variance filtering (0 for uniform variance)
  param(varianceDef.var_niter, varianceConf);
  // Number of passes for the variance filtering (0 for uniform variance)
  param(varianceDef.var_npass, varianceConf);

  // Optimality test section
  OptimalityTestDef optimalityTestDef;
  eckit::LocalConfiguration optimalityTestConf;
  // Number of length-scale factors for optimization
  param(optimalityTestDef.optimality_nfac, optimalityTestConf);
  // Increments of length-scale factors for optimization
  param(optimalityTestDef.optimality_delta, optimalityTestConf);
  // Number of test vectors for optimization
  param(optimalityTestDef.optimality_ntest, optimalityTestConf);

  // Fit section
  FitDef fitDef;
  eckit::LocalConfiguration fitConf;
  // Horizontal filtering suport radius [in meters]
  param(fitDef.diag_rhflt, fitConf);
  // Vertical filtering support radius
  param(fitDef.diag_rvflt, fitConf);
  // Number of levels between interpolation levels
  param(fitDef.fit_dl0, fitConf);
  // Number of components in the fit function
  param(fitDef.fit_ncmp, fitConf);

  // NICAS section
  NICASDef nicasDef;
  eckit::LocalConfiguration nicasConf;
  // Resolution
  param(nicasDef.resol, nicasConf);
  // Maximum size of the Sc1 subset
  param(nicasDef.nc1max, nicasConf);
  // NICAS draw type ('random' or 'octahedral')
  param(nicasDef.nicas_draw_type, nicasConf);
  // Force specific support radii
  param(nicasDef.forced_radii, nicasConf);
  // Normalization randomization size
  param(nicasDef.norm_rand_size, nicasConf);
  // Positive-definiteness test
  param(nicasDef.pos_def_test, nicasConf);
  // Horizontal NICAS interpolation test
  param(nicasDef.interp_test, nicasConf);
  // Overriding component in file
  param(nicasDef.file_component, nicasConf);

  // Psichitouv section
  PsichitouvDef psichitouvDef;
  eckit::LocalConfiguration psichitouvConf;
  // Number of longitudes for the regular grid
  param(psichitouvDef.wind_nlon, psichitouvConf);
  // Number of latitudes for the regular grid
  param(psichitouvDef.wind_nlat, psichitouvConf);
  // Half-width of the Savitzky-Golay to compute derivatives
  param(psichitouvDef.wind_nsg, psichitouvConf);
  // Wind inflation to compensate the Savitzky-Golay smoothing
  param(psichitouvDef.wind_inflation, psichitouvConf);

  // External section
  ExternalDef externalDef;
  eckit::LocalConfiguration externalConf;
  // Iterative algorithm (ensemble members loaded sequentially)
  param(externalDef.iterative_algo, externalConf);

  // General parameters
  config->set("general", generalConf);
  // IO parameters
  config->set("io", ioConf);
  // Drivers parameters
  config->set("drivers", driversConf);
  // Model parameters
  config->set("model", modelConf);
  // Ensemble sizes parameters
  config->set("ensemble sizes", ensembleSizesConf);
  // Sampling parameters
  config->set("sampling", samplingConf);
  // Diagnostics parameters
  config->set("diagnostics", diagnosticsConf);
  // Vertical balance parameters
  config->set("vertical balance", verticalBalanceConf);
  // Variance parameters
  config->set("variance", varianceConf);
  // Optimality test parameters
  config->set("optimality test", optimalityTestConf);
  // Fit parameters
  config->set("fit", fitConf);
  // NICAS parameters
  config->set("nicas", nicasConf);
  // Psichitouv parameters
  config->set("psichitouv", psichitouvConf);
  // External parameters
  config->set("external", externalConf);
}

// -----------------------------------------------------------------------------

}  // namespace bump_lib
