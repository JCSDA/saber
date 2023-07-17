/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <math.h>

#include <limits>
#include <string>
#include <utility>

#include "eckit/config/Configuration.h"

namespace bump_lib {

// -----------------------------------------------------------------------------

// General section
struct GeneralDef {
  // Add colors to the log (for display on terminal)
  std::pair<const char *, bool> color_log =
    std::make_pair("color log", false);

  // Stream test messages into a dedicated channel
  std::pair<const char *, bool> testing =
    std::make_pair("testing", false);

  // Default seed for random numbers
  std::pair<const char *, bool> default_seed =
    std::make_pair("default seed", true);

  // Inter-compilers reproducibility
  std::pair<const char *, bool> repro_ops =
    std::make_pair("reproducibility operators", true);

  // Reproducibility threshold
  std::pair<const char *, double> repro_th =
    std::make_pair("reproducibility threshold", 1.0e-12);

  // Universe radius [in meters]
  std::pair<const char *, double> universe_radius =
    std::make_pair("universe length-scale", 6371229*M_PI);

  // Sampling method
  std::pair<const char *, std::string> sampling_method =
    std::make_pair("sampling method", "potential");
};

// I/O section
struct IODef {
  // Data directory
  std::pair<const char *, std::string> data_directory =
    std::make_pair("data directory", ".");

  // Data prefix
  std::pair<const char *, std::string> files_prefix =
    std::make_pair("files prefix", "");

  // Write in new files
  std::pair<const char *, bool> new_files =
    std::make_pair("new files", true);

  // Parallel NetCDF I/O
  std::pair<const char *, bool> parallel_netcdf =
    std::make_pair("parallel netcdf", true);

  // Number of I/O processors
  std::pair<const char *, int> nprocio =
    std::make_pair("io tasks", 20);

  // Sampling file
  std::pair<const char *, std::string> fname_samp =
    std::make_pair("overriding sampling file", "");

  // Vertical balance file
  std::pair<const char *, std::string> fname_vbal =
    std::make_pair("overriding vertical balance file", "");

  // Universe radius file
  std::pair<const char *, std::string> fname_universe_radius =
    std::make_pair("overriding universe radius file", "");

  // NICAS file
  std::pair<const char *, std::string> fname_nicas =
    std::make_pair("overriding nicas file", "");

  // Psichitouv transform file
  std::pair<const char *, std::string> fname_wind =
    std::make_pair("overriding psichitouv file", "");

  // GSI data file
  std::pair<const char *, std::string> fname_gsi_data =
    std::make_pair("gsi data file", "");

  // GSI namelist
  std::pair<const char *, std::string> fname_gsi_nam =
    std::make_pair("gsi namelist", "");
};

// Drivers section
struct DriversDef {
  // Compute covariance, ensemble 1
  std::pair<const char *, bool> compute_cov1 =
    std::make_pair("compute covariance", false);

  // Compute covariance, ensemble 2
  std::pair<const char *, bool> compute_cov2 =
    std::make_pair("compute lowres covariance", false);

  // Compute correlation, ensemble 1
  std::pair<const char *, bool> compute_cor1 =
    std::make_pair("compute correlation", false);

  // Compute correlation, ensemble 2
  std::pair<const char *, bool> compute_cor2 =
    std::make_pair("compute lowres correlation", false);

  // Compute localization, ensemble 1
  std::pair<const char *, bool> compute_loc1 =
    std::make_pair("compute localization", false);

  // Compute localization, ensemble 2
  std::pair<const char *, bool> compute_loc2 =
    std::make_pair("compute lowres localization", false);

  // Compute hybrid weights
  std::pair<const char *, bool> compute_hyb =
    std::make_pair("compute hybrid weights", false);

  // Hybrid term source ('randomized static' or 'lowres ensemble')
  std::pair<const char *, std::string> hybrid_source =
    std::make_pair("hybrid source", "");

  // Multivariate strategy ('univariate', 'duplicated', 'duplicated and weighted' or 'crossed')
  std::pair<const char *, std::string> strategy =
    std::make_pair("multivariate strategy", "");

  // New normality test
  std::pair<const char *, bool> new_normality =
    std::make_pair("compute normality", false);

  // Read local sampling
  std::pair<const char *, bool> load_samp_local =
    std::make_pair("read local sampling", false);

  // Read global sampling
  std::pair<const char *, bool> load_samp_global =
    std::make_pair("read global sampling", false);

  // Write local sampling
  std::pair<const char *, bool> write_samp_local =
    std::make_pair("write local sampling", false);

  // Write global sampling
  std::pair<const char *, bool> write_samp_global =
    std::make_pair("write global sampling", false);

  // Write sampling grids
  std::pair<const char *, bool> write_samp_grids =
    std::make_pair("write sampling grids", false);

  // New vertical covariance
  std::pair<const char *, bool> new_vbal_cov =
    std::make_pair("compute vertical covariance", false);

  // Read local vertical covariance
  std::pair<const char *, bool> load_vbal_cov =
    std::make_pair("read vertical covariance", false);

  // Write local vertical covariancee
  std::pair<const char *, bool> write_vbal_cov =
    std::make_pair("write vertical covariance", false);

  // Compute vertical balance operator
  std::pair<const char *, bool> new_vbal =
    std::make_pair("compute vertical balance", false);

  // Read local vertical balance operator
  std::pair<const char *, bool> load_vbal =
    std::make_pair("read vertical balance", false);

  // Write vertical balance operator
  std::pair<const char *, bool> write_vbal =
    std::make_pair("write vertical balance", false);

  // Compute variance
  std::pair<const char *, bool> new_var =
    std::make_pair("compute variance", false);

  // Compute moments
  std::pair<const char *, bool> new_mom =
    std::make_pair("compute moments", false);

  // Read sampling moments
  std::pair<const char *, bool> load_mom =
    std::make_pair("read moments", false);

  // Write sampling moments
  std::pair<const char *, bool> write_mom =
    std::make_pair("write moments", false);

  // Write HDIAG diagnostics
  std::pair<const char *, bool> write_hdiag =
    std::make_pair("write diagnostics", false);

  // Write HDIAG components detail
  std::pair<const char *, bool> write_hdiag_detail =
    std::make_pair("write diagnostics detail", false);

  // Read universe radius
  std::pair<const char *, bool> load_universe_radius =
    std::make_pair("read universe radius", false);

  // Write universe radius
  std::pair<const char *, bool> write_universe_radius =
    std::make_pair("write universe radius", false);

  // Compute NICAS
  std::pair<const char *, bool> new_nicas =
    std::make_pair("compute nicas", false);

  // Read local NICAS parameters
  std::pair<const char *, bool> load_nicas_local =
    std::make_pair("read local nicas", false);

  // Read global NICAS parameters
  std::pair<const char *, bool> load_nicas_global =
    std::make_pair("read global nicas", false);

  // Write local NICAS parameters
  std::pair<const char *, bool> write_nicas_local =
    std::make_pair("write local nicas", false);

  // Write global NICAS parameters
  std::pair<const char *, bool> write_nicas_global =
    std::make_pair("write global nicas", false);

  // Write NICAS grids
  std::pair<const char *, bool> write_nicas_grids =
    std::make_pair("write nicas grids", false);

  // Write NICAS steps
  std::pair<const char *, bool> write_nicas_steps =
    std::make_pair("write nicas steps", false);

  // Compute wind transform
  std::pair<const char *, bool> new_wind =
    std::make_pair("compute psichitouv", false);

  // Read local wind transform
  std::pair<const char *, bool> load_wind_local =
    std::make_pair("read local psichitouv", false);

  // Write local wind transform
  std::pair<const char *, bool> write_wind_local =
    std::make_pair("write local psichitouv", false);

  // Test vertical balance inverse
  std::pair<const char *, bool> check_vbal =
    std::make_pair("vertical balance inverse test", false);

  // Test adjoints
  std::pair<const char *, bool> check_adjoints =
    std::make_pair("adjoints test", false);

  // Test NICAS normalization (number of tests)
  std::pair<const char *, int> check_normalization =
    std::make_pair("normalization test", 0);

  // Test NICAS application on diracs
  std::pair<const char *, bool> check_dirac =
    std::make_pair("internal dirac test", false);

  // Test NICAS randomization
  std::pair<const char *, bool> check_randomization =
    std::make_pair("randomization test", false);

  // Test HDIAG-NICAS consistency
  std::pair<const char *, bool> check_consistency =
    std::make_pair("internal consistency test", false);

  // Test HDIAG optimality
  std::pair<const char *, bool> check_optimality =
    std::make_pair("localization optimality test", false);

  // Interpolate vertical balance, standard-deviation or length-scales from GSI data
  std::pair<const char *, bool> from_gsi =
    std::make_pair("interpolate from gsi data", false);
};

// Model section
struct ModelDef {
  // Level for 2D variables ('first' or 'last')
  std::pair<const char *, std::string> lev2d =
    std::make_pair("level for 2d variables", "first");

  // Check that sampling couples and interpolations do not cross mask boundaries
  std::pair<const char *, bool> mask_check =
    std::make_pair("do not cross mask boundaries", false);
};

// Ensemble sizes section
struct EnsembleSizesDef {
  // Ensemble 1 size
  std::pair<const char *, int> ens1_ne =
    std::make_pair("total ensemble size", 0);

  // Ensemble 1 sub-ensembles number
  std::pair<const char *, int> ens1_nsub =
    std::make_pair("sub-ensembles", 1);

  // Ensemble 2 size
  std::pair<const char *, int> ens2_ne =
    std::make_pair("total lowres ensemble size", 0);

  // Ensemble 2 sub-ensembles number
  std::pair<const char *, int> ens2_nsub =
    std::make_pair("lowres sub-ensembles", 1);
};

// Mask parameters
struct MaskDef {
  // Mask threshold
  std::pair<const char *, double> mask_th =
    std::make_pair("threshold", 0.0);

  // Mask threshold side ('lower' if mask_th is the lower bound, resp. 'upper')
  std::pair<const char *, std::string> mask_lu =
    std::make_pair("side", "");

  // Mask variable
  std::pair<const char *, std::string> mask_variable =
    std::make_pair("variable", "");
};

// Sampling section
struct SamplingDef {
  // Computation grid size
  std::pair<const char *, int> nc1 =
    std::make_pair("computation grid size", 0);

  // Diagnostic grid size
  std::pair<const char *, int> nc2 =
    std::make_pair("diagnostic grid size", 0);

  // Number of distance classes
  std::pair<const char *, int> nc3 =
    std::make_pair("distance classes", 0);

  // Number of angular sectors
  std::pair<const char *, int> nc4 =
    std::make_pair("angular sectors", 1);

  // Class size (for sam_type='hor'), should be larger than the typical grid cell size [in meters]
  std::pair<const char *, double> dc =
    std::make_pair("distance class width", 0.0);

  // Reduced number of levels for diagnostics
  std::pair<const char *, int> nl0r =
    std::make_pair("reduced levels", 0);

  // Activate local diagnostics
  std::pair<const char *, bool> local_diag =
    std::make_pair("local diagnostic", false);

  // Local diagnostics calculation radius [in meters]
  std::pair<const char *, double> local_rad =
    std::make_pair("averaging length-scale", 0.0);

  // Local diagnostics calculation latitude band half-width [in degrees]
  std::pair<const char *, double> local_dlat =
    std::make_pair("averaging latitude width", 0.0);

  // Diagnostic draw type ('random' or 'octahedral')
  std::pair<const char *, std::string> draw_type =
    std::make_pair("grid type", "random");

  // Maximum number of random number draws
  std::pair<const char *, int> irmax =
    std::make_pair("max number of draws", 10000);

  // Vertical balance C2B to C0A interpolation type ('c0': C0 mesh-based, 'c1': C1 mesh-based
  // or 'si': smooth interpolation)
  std::pair<const char *, std::string> interp_type =
    std::make_pair("interpolation type", "si");

  // Threshold on vertically contiguous points for the mask (0 to skip the test)
  std::pair<const char *, int> ncontig_th =
    std::make_pair("contiguous levels threshold", 0);
};

// Diagnostics section
struct DiagnosticsDef {
  // Ensemble size
  std::pair<const char *, int> ne =
    std::make_pair("target ensemble size", 0);

  // Ensemble size of the hybrid term
  std::pair<const char *, int> ne_lr =
    std::make_pair("target lowres ensemble size", 0);

  // Gaussian approximation for asymptotic quantities
  std::pair<const char *, bool> gau_approx =
    std::make_pair("gaussian approximation", false);

  // Threshold on generalized kurtosis (3.0 = Gaussian distribution)
  std::pair<const char *, double> gen_kurt_th =
    std::make_pair("generalized kurtosis threshold", std::numeric_limits<double>().max());

  // Number of bins for averaged statistics histograms
  std::pair<const char *, int> avg_nbins =
    std::make_pair("histogram bins", 0);

  // Support radius scaling in CMAT from HDIAG
  std::pair<const char *, double> lengths_scaling =
    std::make_pair("diagnosed lengths scaling", 1.0);
};

// Vertical balance block parameters
struct VerticalBalanceBlockDef {
  // Diagonal auto-covariance for the inversion
  std::pair<const char *, bool> diag_auto =
    std::make_pair("diagonal autocovariance", false);

  // Diagonal regression
  std::pair<const char *, bool> diag_reg =
    std::make_pair("diagonal regression", false);

  // Scalar coefficients for identity vertical balance
  std::pair<const char *, double> id_coef =
    std::make_pair("identity block weight", 1.0);
};

// Vertical balance section
struct VerticalBalanceDef {
  // Pseudo-inverse for auto-covariance
  std::pair<const char *, bool> vbal_pseudo_inv =
    std::make_pair("pseudo inverse", false);

  // Dominant mode for pseudo-inverse
  std::pair<const char *, int> vbal_pseudo_inv_mmax =
    std::make_pair("dominant mode", 0);

  // Variance threshold to compute the dominant mode for pseudo-inverse
  std::pair<const char *, double> vbal_pseudo_inv_var_th =
    std::make_pair("variance threshold", 0.0);

  // Identity vertical balance for tests
  std::pair<const char *, bool> vbal_id =
    std::make_pair("identity blocks", false);
};

// Variance section
struct VarianceDef {
  // Force specific variance
  std::pair<const char *, bool> forced_var =
    std::make_pair("explicit stddev", false);

  // Filter variance
  std::pair<const char *, bool> var_filter =
    std::make_pair("objective filtering", false);

  // Number of iterations for the variance filtering (0 for uniform variance)
  std::pair<const char *, int> var_niter =
    std::make_pair("filtering iterations", -1);

  // Number of passes for the variance filtering (0 for uniform variance)
  std::pair<const char *, int> var_npass =
    std::make_pair("filtering passes", -1);
};

// Optimality test section
struct OptimalityTestDef {
  // Number of length-scale factors for optimization
  std::pair<const char *, int> optimality_nfac =
    std::make_pair("half number of factors", 1);

  // Increments of length-scale factors for optimization
  std::pair<const char *, double> optimality_delta =
    std::make_pair("factors increment", 0.05);

  // Number of test vectors for optimization
  std::pair<const char *, int> optimality_ntest =
    std::make_pair("test vectors", 10);
};

// Fit section
struct FitDef {
  // Horizontal filtering suport radius [in meters]
  std::pair<const char *, double> diag_rhflt =
    std::make_pair("horizontal filtering length-scale", 0.0);

  // Vertical filtering support radius
  std::pair<const char *, double> diag_rvflt =
    std::make_pair("vertical filtering length-scale", 0.0);

  // Number of levels between interpolation levels
  std::pair<const char *, int> fit_dl0 =
    std::make_pair("vertical stride", 1);

  // Number of components in the fit function
  std::pair<const char *, int> fit_ncmp =
    std::make_pair("number of components", 1);
};

// NICAS section
struct NICASDef {
  // Resolution
  std::pair<const char *, double> resol =
    std::make_pair("resolution", 0.0);

  // Maximum size of the Sc1 subset
  std::pair<const char *, int> nc1max =
    std::make_pair("max horizontal grid size", 15000);

  // NICAS draw type ('random' or 'octahedral')
  std::pair<const char *, std::string> nicas_draw_type =
    std::make_pair("grid type", "random");

  // Force specific support radii
  std::pair<const char *, bool> forced_radii =
    std::make_pair("explicit length-scales", false);

  // Normalization randomization size
  std::pair<const char *, int> norm_rand_size =
    std::make_pair("normalization randomization size", 0);

  // Positive-definiteness test
  std::pair<const char *, bool> pos_def_test =
    std::make_pair("positive-definiteness test", false);

  // Horizontal NICAS interpolation test
  std::pair<const char *, bool> interp_test =
    std::make_pair("horizontal interpolation test", false);

  // Overriding component in file
  std::pair<const char *, int> file_component =
    std::make_pair("overriding component in file", 0);
};

// Psichitouv section
struct PsichitouvDef {
  // Number of longitudes for the regular grid
  std::pair<const char *, int> wind_nlon =
    std::make_pair("longitudes", 0);

  // Number of latitudes for the regular grid
  std::pair<const char *, int> wind_nlat =
    std::make_pair("latitudes", 0);

  // Half-width of the Savitzky-Golay to compute derivatives
  std::pair<const char *, int> wind_nsg =
    std::make_pair("savitzky-golay half width", 0);

  // Wind inflation to compensate the Savitzky-Golay smoothing
  std::pair<const char *, double> wind_inflation =
    std::make_pair("wind inflation", 1.0);
};

// External section
struct ExternalDef {
  // Iterative algorithm (ensemble members loaded sequentially)
  std::pair<const char *, bool> iterative_algo =
    std::make_pair("iterative algorithm", false);
};

// -----------------------------------------------------------------------------

extern "C" {
  void bump_config_init_f90(eckit::LocalConfiguration *);
}

// -----------------------------------------------------------------------------

}  // namespace bump_lib
