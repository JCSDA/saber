/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/util/missingValues.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/bump/type_bump_parameters.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------
// Elemental parameters (without default value)
// -----------------------------------------------------------------------------

// Variables value or profile elemental parameters
class VarsValueOrProfileParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VarsValueOrProfileParameters, oops::Parameters)

 public:
  // Variables
  oops::RequiredParameter<std::vector<std::string>> variables{"variables", this};
  // Value
  oops::OptionalParameter<double> value{"value", this};
  // Profile
  oops::OptionalParameter<std::vector<double>> profile{"profile", this};
};

// -----------------------------------------------------------------------------

// Groups value or profile elemental parameters
class GroupsValueOrProfileParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GroupsValueOrProfileParameters, oops::Parameters)

 public:
  // Groups
  oops::RequiredParameter<std::vector<std::string>> groups{"groups", this};
  // Value
  oops::OptionalParameter<double> value{"value", this};
  // Profile
  oops::OptionalParameter<std::vector<double>> profile{"profile", this};
};

// -----------------------------------------------------------------------------

// Groups value elemental parameters
class GroupsValueParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GroupsValueParameters, oops::Parameters)

 public:
  // Variables
  oops::RequiredParameter<std::vector<std::string>> groups{"groups", this};
  // Value
  oops::RequiredParameter<int> value{"value", this};
};

// -----------------------------------------------------------------------------

// Alias elemental paramaters
class AliasParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(AliasParameters, oops::Parameters)

 public:
  // In code
  oops::RequiredParameter<std::string> in_code{"in code", this};
  // In file
  oops::RequiredParameter<std::string> in_file{"in file", this};
};

// -----------------------------------------------------------------------------

class GroupParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GroupParameters, oops::Parameters)

 public:
  // Group name
  oops::RequiredParameter<std::string> name{"group name", this};
  // Group variables
  oops::RequiredParameter<std::vector<std::string>> variables{"variables", this};
};

// -----------------------------------------------------------------------------

// Local profile elemental parameters
class LocalProfileParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LocalProfileParameters, oops::Parameters)

 public:
  // Longitudes of the local diagnostics profiles to write [in degrees]
  oops::RequiredParameter<double> lon_ldwv{"longitude", this};
  // Latitudes of the local diagnostics profiles to write [in degrees]
  oops::RequiredParameter<double> lat_ldwv{"latitude", this};
  // Name of the local diagnostics profiles to write
  oops::RequiredParameter<std::string> name_ldwv{"name", this};
};

// -----------------------------------------------------------------------------

// Groups type elemental parameters
class GroupsTypeParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GroupsTypeParameters, oops::Parameters)

 public:
  // Groups
  oops::RequiredParameter<std::vector<std::string>> groups{"groups", this};
  // Type
  oops::RequiredParameter<std::string> type{"type", this};
};

// -----------------------------------------------------------------------------

// Local weight elemental parameters
class LocWgtParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LocWgtParameters, oops::Parameters)

 public:
  // Row variables
  oops::RequiredParameter<std::vector<std::string>> row_variables{"row variables", this};
  // Column variables
  oops::RequiredParameter<std::vector<std::string>> column_variables{"column variables", this};
  // Value
  oops::RequiredParameter<double> value{"value", this};
};

// -----------------------------------------------------------------------------

// Dirac point elemental parameters
class DiracPointParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DiracPointParameters, oops::Parameters)

 public:
  // Diracs longitudes [in degrees]
  oops::RequiredParameter<double> longitude{"longitude", this};
  // Diracs latitudes [in degrees]
  oops::RequiredParameter<double> latitude{"latitude", this};
  // Diracs level
  oops::RequiredParameter<int> level{"level", this};
  // Diracs variable indices
  oops::RequiredParameter<std::string> variable{"variable", this};
};

// -----------------------------------------------------------------------------
// Utility
// -----------------------------------------------------------------------------

// Shorcut to create a Parameter based on a pair <name, default value>
template <class T>
oops::Parameter<T> param(std::pair<const char*, T> pair, oops::Parameters *parent) {
  oops::Parameter<T> par{pair.first, pair.second, parent};
  return par;
}

// -----------------------------------------------------------------------------
// BUMP parameters sections
// -----------------------------------------------------------------------------

// General section
class GeneralSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeneralSection, oops::Parameters)

 private:
  GeneralDef def;

 public:
  // Add colors to the log (for display on terminal)
  oops::Parameter<bool> color_log = param(def.color_log, this);
  // Stream test messages into a dedicated channel
  oops::Parameter<bool> testing = param(def.testing, this);
  // Default seed for random numbers (0 for time-dependent seed)
  oops::Parameter<int> default_seed = param(def.default_seed, this);
  // Inter-compilers reproducibility
  oops::Parameter<bool> repro_ops = param(def.repro_ops, this);
  // Reproducibility threshold
  oops::Parameter<double> repro_th = param(def.repro_th, this);
  // Universe radius [in meters]
  oops::Parameter<double> universe_radius = param(def.universe_radius, this);
  // Sampling method
  oops::Parameter<std::string> sampling_method = param(def.sampling_method, this);
};

// -----------------------------------------------------------------------------

// I/O section
class IOSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(IOSection, oops::Parameters)

 private:
  IODef def;

 public:
  // Data directory
  oops::Parameter<std::string> data_directory = param(def.data_directory, this);
  // Data prefix
  oops::Parameter<std::string> files_prefix = param(def.files_prefix, this);
  // Write in new files
  oops::Parameter<bool> new_files = param(def.new_files, this);
  // Parallel NetCDF I/O
  oops::Parameter<bool> parallel_netcdf = param(def.parallel_netcdf, this);
  // Number of I/O processors
  oops::Parameter<int> nprocio = param(def.nprocio, this);
  // Alias
  oops::Parameter<std::vector<AliasParameters>> alias{"alias", {}, this};
  // Sampling file
  oops::Parameter<std::string> fname_samp = param(def.fname_samp, this);
  // Vertical covariance files
  oops::Parameter<std::vector<std::string>> fname_vbal_cov{"overriding vertical covariance file",
    {}, this};
  // Vertical balance file
  oops::Parameter<std::string> fname_vbal = param(def.fname_vbal, this);
  // Ensemble 1 moments files
  oops::Parameter<std::vector<std::string>> fname_mom{"overriding moments file", {}, this};
  // Ensemble 2 moments files
  oops::Parameter<std::vector<std::string>> fname_mom2{"overriding lowres moments file", {}, this};
  // Universe radius file
  oops::Parameter<std::string> fname_universe_radius = param(def.fname_universe_radius, this);
  // NICAS file
  oops::Parameter<std::string> fname_nicas = param(def.fname_nicas, this);
  // Psichitouv transform file
  oops::Parameter<std::string> fname_wind = param(def.fname_wind, this);
  // GSI data file
  oops::Parameter<std::string> fname_gsi_data = param(def.fname_gsi_data, this);
  // GSI namelist
  oops::Parameter<std::string> fname_gsi_nam = param(def.fname_gsi_nam, this);
};

// -----------------------------------------------------------------------------

// Drivers section
class DriversSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DriversSection, oops::Parameters)

 private:
  DriversDef def;

 public:
  // Compute covariance, ensemble 1
  oops::Parameter<bool> compute_cov1 = param(def.compute_cov1, this);
  // Compute covariance, ensemble 2
  oops::Parameter<bool> compute_cov2 = param(def.compute_cov2, this);
  // Compute correlation, ensemble 1
  oops::Parameter<bool> compute_cor1 = param(def.compute_cor1, this);
  // Compute correlation, ensemble 2
  oops::Parameter<bool> compute_cor2 = param(def.compute_cor2, this);
  // Compute localization, ensemble 1
  oops::Parameter<bool> compute_loc1 = param(def.compute_loc1, this);
  // Compute localization, ensemble 2
  oops::Parameter<bool> compute_loc2 = param(def.compute_loc2, this);
  // Compute hybrid weights
  oops::Parameter<bool> compute_hyb = param(def.compute_hyb, this);
  // Hybrid term source ('randomized static' or 'lowres ensemble')
  oops::Parameter<std::string> hybrid_source = param(def.hybrid_source, this);
  // Multivariate strategy ('univariate', 'duplicated', 'duplicated and weighted' or 'crossed')
  oops::Parameter<std::string> strategy = param(def.strategy, this);
  // New normality test
  oops::Parameter<bool> new_normality = param(def.new_normality, this);
  // Read local sampling
  oops::Parameter<bool> load_samp_local = param(def.load_samp_local, this);
  // Read global sampling
  oops::Parameter<bool> load_samp_global = param(def.load_samp_global, this);
  // Write local sampling
  oops::Parameter<bool> write_samp_local = param(def.write_samp_local, this);
  // Write global sampling
  oops::Parameter<bool> write_samp_global = param(def.write_samp_global, this);
  // Write sampling grids
  oops::Parameter<bool> write_samp_grids = param(def.write_samp_grids, this);
  // New vertical covariance
  oops::Parameter<bool> new_vbal_cov = param(def.new_vbal_cov, this);
  // Read local vertical covariance
  oops::Parameter<bool> load_vbal_cov = param(def.load_vbal_cov, this);
  // Write local vertical covariance
  oops::Parameter<bool> write_vbal_cov = param(def.write_vbal_cov, this);
  // Compute vertical balance operator
  oops::Parameter<bool> new_vbal = param(def.new_vbal, this);
  // Read local vertical balance operator
  oops::Parameter<bool> load_vbal = param(def.load_vbal, this);
  // Write vertical balance operator
  oops::Parameter<bool> write_vbal = param(def.write_vbal, this);
  // Compute variance
  oops::Parameter<bool> new_var = param(def.new_var, this);
  // Compute moments
  oops::Parameter<bool> new_mom = param(def.new_mom, this);
  // Read sampling moments
  oops::Parameter<bool> load_mom = param(def.load_mom, this);
  // Write sampling moments
  oops::Parameter<bool> write_mom = param(def.write_mom, this);
  // Write HDIAG diagnostics
  oops::Parameter<bool> write_hdiag = param(def.write_hdiag, this);
  // Write HDIAG components detail
  oops::Parameter<bool> write_hdiag_detail = param(def.write_hdiag_detail, this);
  // Read universe radius
  oops::Parameter<bool> load_universe_radius = param(def.load_universe_radius, this);
  // Write universe radius
  oops::Parameter<bool> write_universe_radius = param(def.write_universe_radius, this);
  // Compute NICAS
  oops::Parameter<bool> new_nicas = param(def.new_nicas, this);
  // Read local NICAS parameters
  oops::Parameter<bool> load_nicas_local = param(def.load_nicas_local, this);
  // Read global NICAS parameters
  oops::Parameter<bool> load_nicas_global = param(def.load_nicas_global, this);
  // Write local NICAS parameters
  oops::Parameter<bool> write_nicas_local = param(def.write_nicas_local, this);
  // Write global NICAS parameters
  oops::Parameter<bool> write_nicas_global = param(def.write_nicas_global, this);
  // Write NICAS grids
  oops::Parameter<bool> write_nicas_grids = param(def.write_nicas_grids, this);
  // Write NICAS steps
  oops::Parameter<bool> write_nicas_steps = param(def.write_nicas_steps, this);
  // Compute wind transform
  oops::Parameter<bool> new_wind = param(def.new_wind, this);
  // Read local wind transform
  oops::Parameter<bool> load_wind_local = param(def.load_wind_local, this);
  // Write local wind transform
  oops::Parameter<bool> write_wind_local = param(def.write_wind_local, this);
  // Test vertical balance inverse
  oops::Parameter<bool> check_vbal = param(def.check_vbal, this);
  // Test adjoints
  oops::Parameter<bool> check_adjoints = param(def.check_adjoints, this);
  // Test NICAS normalization (number of tests)
  oops::Parameter<int> check_normalization = param(def.check_normalization, this);
  // Test NICAS application on diracs
  oops::Parameter<bool> check_dirac = param(def.check_dirac, this);
  // Test NICAS randomization
  oops::Parameter<bool> check_randomization = param(def.check_randomization, this);
  // Test HDIAG-NICAS consistency
  oops::Parameter<bool> check_consistency = param(def.check_consistency, this);
  // Test HDIAG optimality
  oops::Parameter<bool> check_optimality = param(def.check_optimality, this);
  // Interpolate vertical balance, standard-deviation or length-scales from GSI data
  oops::Parameter<bool> from_gsi = param(def.from_gsi, this);
};

// -----------------------------------------------------------------------------

// Model section
class ModelSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelSection, oops::Parameters)

 private:
  ModelDef def;

 public:
  // Level for 2D variables ('first' or 'last')
  oops::Parameter<std::string> lev2d = param(def.lev2d, this);
  // Variables names
  oops::Parameter<std::vector<std::string>> variables{"variables", {}, this};
  // 2D variables names
  oops::Parameter<std::vector<std::string>> var2d{"2d variables", {}, this};
  // Groups of variables
  oops::OptionalParameter<std::vector<GroupParameters>> groups{"groups", this};
  // Check that sampling couples and interpolations do not cross mask boundaries
  oops::Parameter<bool> mask_check = param(def.mask_check, this);
};

// -----------------------------------------------------------------------------

// Ensemble size section
class EnsembleSizesSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(EnsembleSizesSection, oops::Parameters)

 private:
  EnsembleSizesDef def;

 public:
  // Ensemble 1 size
  oops::Parameter<int> ens1_ne = param(def.ens1_ne, this);
  // Ensemble 1 sub-ensembles number
  oops::Parameter<int> ens1_nsub = param(def.ens1_nsub, this);
  // Ensemble 2 size
  oops::Parameter<int> ens2_ne = param(def.ens2_ne, this);
  // Ensemble 2 sub-ensembles number
  oops::Parameter<int> ens2_nsub = param(def.ens2_nsub, this);
};

// -----------------------------------------------------------------------------

// Mask parameters
class MaskParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(MaskParameters, oops::Parameters)

 private:
  MaskDef def;

 public:
  // Mask restriction type
  oops::RequiredParameter<std::string> mask_type{"type", this};
  // Mask threshold
  oops::Parameter<double> mask_th = param(def.mask_th, this);
  // Mask threshold side ('lower' if mask_th is the lower bound, resp. 'upper')
  oops::Parameter<std::string> mask_lu = param(def.mask_lu, this);
  // Mask variable
  oops::Parameter<std::string> mask_variable = param(def.mask_variable, this);
};

// Sampling section
class SamplingSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SamplingSection, oops::Parameters)

 private:
  SamplingDef def;

 public:
  // Computation grid size
  oops::Parameter<int> nc1 = param(def.nc1, this);
  // Diagnostic grid size
  oops::Parameter<int> nc2 = param(def.nc2, this);
  // Number of distance classes
  oops::Parameter<int> nc3 = param(def.nc3, this);
  // Number of angular sectors
  oops::Parameter<int> nc4 = param(def.nc4, this);
  // Class size (for sam_type='hor'), should be larger than the typical grid cell size [in meters]
  oops::Parameter<double> dc = param(def.dc, this);
  // Reduced number of levels for diagnostics
  oops::Parameter<int> nl0r = param(def.nl0r, this);
  // Activate local diagnostics
  oops::Parameter<bool> local_diag = param(def.local_diag, this);
  // Local diagnostics calculation radius [in meters]
  oops::Parameter<double> local_rad = param(def.local_rad, this);
  // Local diagnostics calculation latitude band half-width [in degrees]
  oops::Parameter<double> local_dlat = param(def.local_dlat, this);
  // Diagnostic draw type ('random' or 'octahedral')
  oops::Parameter<std::string> draw_type = param(def.draw_type, this);
  // Maximum number of random number draws
  oops::Parameter<int> irmax = param(def.irmax, this);
  // Vertical balance C2B to C0A interpolation type ('c0': C0 mesh-based, 'c1': C1 mesh-based
  // or 'si': smooth interpolation)
  oops::Parameter<std::string> interp_type = param(def.interp_type, this);
  // Sampling masks
  oops::Parameter<std::vector<MaskParameters>> masks{"masks", {}, this};
  // Threshold on vertically contiguous points for the mask (0 to skip the test)
  oops::Parameter<int> ncontig_th = param(def.ncontig_th, this);
};

// -----------------------------------------------------------------------------

// Diagnostics section
class DiagnosticsSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DiagnosticsSection, oops::Parameters)

 private:
  DiagnosticsDef def;

 public:
  // Ensemble size
  oops::Parameter<int> ne = param(def.ne, this);
  // Ensemble size of the hybrid term
  oops::Parameter<int> ne_lr = param(def.ne_lr, this);
  // Gaussian approximation for asymptotic quantities
  oops::Parameter<bool> gau_approx = param(def.gau_approx, this);
  // Compute localization from correlation
  oops::Parameter<bool> loc_from_cor = param(def.loc_from_cor, this);
  // Threshold on generalized kurtosis (3.0 = Gaussian distribution)
  oops::Parameter<double> gen_kurt_th = param(def.gen_kurt_th, this);
  // Number of bins for averaged statistics histograms
  oops::Parameter<int> avg_nbins = param(def.avg_nbins, this);
  // Support radius scaling in CMAT from HDIAG
  oops::Parameter<double> lengths_scaling = param(def.lengths_scaling, this);
};

// -----------------------------------------------------------------------------

// Vertical balance block parameters
class VerticalBalanceBlockParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VerticalBalanceBlockParameters, oops::Parameters)

 private:
  VerticalBalanceBlockDef def;

 public:
  // Balanced variable
  oops::RequiredParameter<std::string> balanced{"balanced variable", this};
  // Unbalanced variable
  oops::RequiredParameter<std::string> unbalanced{"unbalanced variable", this};
  // Diagonal auto-covariance for the inversion
  oops::Parameter<bool> diag_auto = param(def.diag_auto, this);
  // Diagonal regression
  oops::Parameter<bool> diag_reg = param(def.diag_reg, this);
  // Scalar coefficients for identity vertical balance
  oops::Parameter<double> id_coef = param(def.id_coef, this);
};

// Vertical balance section
class VerticalBalanceSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VerticalBalanceSection, oops::Parameters)

 private:
  VerticalBalanceDef def;

 public:
  // Vertical balance parameters
  oops::Parameter<std::vector<VerticalBalanceBlockParameters>> vbal{"vbal", {}, this};
  // Pseudo-inverse for auto-covariance
  oops::Parameter<bool> vbal_pseudo_inv = param(def.vbal_pseudo_inv, this);
  // Dominant mode for pseudo-inverse
  oops::Parameter<int> vbal_pseudo_inv_mmax = param(def.vbal_pseudo_inv_mmax, this);
  // Variance threshold to compute the dominant mode for pseudo-inverse
  oops::Parameter<double> vbal_pseudo_inv_var_th = param(def.vbal_pseudo_inv_var_th, this);
  // Identity vertical balance for tests
  oops::Parameter<bool> vbal_id = param(def.vbal_id, this);
};

// -----------------------------------------------------------------------------

// Variance section
class VarianceSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VarianceSection, oops::Parameters)

 private:
  VarianceDef def;

 public:
  // Force specific variance
  oops::Parameter<bool> forced_var = param(def.forced_var, this);
  // Forced standard-deviation
  oops::Parameter<std::vector<VarsValueOrProfileParameters>> stddev{"stddev", {}, this};
  // Filter variance
  oops::Parameter<bool> var_filter = param(def.var_filter, this);
  // Number of iterations for the variance filtering (0 for uniform variance)
  oops::Parameter<int> var_niter = param(def.var_niter, this);
  // Number of passes for the variance filtering (0 for uniform variance)
  oops::Parameter<int> var_npass = param(def.var_npass, this);
  // Variance initial filtering support radius [in meters]
  oops::Parameter<std::vector<VarsValueOrProfileParameters>> var_rhflt{"initial length-scale", {},
    this};
  // Resolution for the NICAS smoother
  oops::Parameter<double> smoother_resol = param(def.smoother_resol, this);
  // Maximum size of the Sc1 subset for the NICAS smoother
  oops::Parameter<int> smoother_nc1max = param(def.smoother_nc1max, this);
  // Minimum effective resolution for the NICAS smoother
  oops::Parameter<double> smoother_resol_eff_min = param(def.smoother_resol_eff_min, this);
};

// -----------------------------------------------------------------------------

// Optimality test section
class OptimalityTestSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OptimalityTestSection, oops::Parameters)

 private:
  OptimalityTestDef def;

 public:
  // Number of length-scale factors for optimization
  oops::Parameter<int> optimality_nfac = param(def.optimality_nfac, this);
  // Increments of length-scale factors for optimization
  oops::Parameter<double> optimality_delta = param(def.optimality_delta, this);
  // Number of test vectors for optimization
  oops::Parameter<int> optimality_ntest = param(def.optimality_ntest, this);
};

// -----------------------------------------------------------------------------

// Fit section
class FitSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FitSection, oops::Parameters)

 private:
  FitDef def;

 public:
  // Threshold to filter out lower raw values
  oops::Parameter<double> diag_raw_th = param(def.diag_raw_th, this);
  // Horizontal filtering suport radius [in meters]
  oops::Parameter<double> diag_rhflt = param(def.diag_rhflt, this);
  // Vertical filtering support radius
  oops::Parameter<double> diag_rvflt = param(def.diag_rvflt, this);
  // Number of levels between interpolation levels
  oops::Parameter<int> fit_dl0 = param(def.fit_dl0, this);
  // Number of components in the fit function
  oops::Parameter<int> fit_ncmp = param(def.fit_ncmp, this);
};

// -----------------------------------------------------------------------------

// NICAS section
class NICASSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(NICASSection, oops::Parameters)

 private:
  NICASDef def;

 public:
  // Resolution
  oops::Parameter<double> resol = param(def.resol, this);
  // Maximum size of the Sc1 subset
  oops::Parameter<int> nc1max = param(def.nc1max, this);
  // Minimum effective resolution
  oops::Parameter<double> resol_eff_min = param(def.resol_eff_min, this);
  // NICAS draw type ('random' or 'octahedral')
  oops::Parameter<std::string> nicas_draw_type = param(def.nicas_draw_type, this);
  // Force specific support radii
  oops::Parameter<bool> forced_radii = param(def.forced_radii, this);
  // Forced horizontal support radius [in meters]
  oops::Parameter<std::vector<GroupsValueOrProfileParameters>> rh{"horizontal length-scale", {},
    this};
  // Forced vertical support radius
  oops::Parameter<std::vector<GroupsValueOrProfileParameters>> rv{"vertical length-scale", {},
    this};
  // Forced localization weights
  oops::Parameter<std::vector<LocWgtParameters>> loc_wgt{"common localization weights", {},
    this};
  // NICAS C1B to C0A interpolation type ('c0': C0 mesh-based, 'c1': C1 mesh-based
  // or 'si': smooth interpolation)
  oops::Parameter<std::vector<GroupsTypeParameters>> interp_type{"interpolation type", {},
    this};
  // Normalization randomization size
  oops::Parameter<int> norm_rand_size = param(def.norm_rand_size, this);
  // Positive-definiteness test
  oops::Parameter<bool> pos_def_test = param(def.pos_def_test, this);
  // Horizontal NICAS interpolation test
  oops::Parameter<bool> interp_test = param(def.interp_test, this);
  // Overriding component in file
  oops::Parameter<int> file_component = param(def.file_component, this);
  // Same horizontal convolution for all levels, no vertical convolution
  oops::Parameter<bool> same_horizontal = param(def.same_horizontal, this);
};

// -----------------------------------------------------------------------------

// Psichitouv section
class PsichitouvSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(PsichitouvSection, oops::Parameters)

 private:
  PsichitouvDef def;

 public:
  // Number of longitudes for the regular grid
  oops::Parameter<int> wind_nlon = param(def.wind_nlon, this);
  // Number of latitudes for the regular grid
  oops::Parameter<int> wind_nlat = param(def.wind_nlat, this);
  // Half-width of the Savitzky-Golay to compute derivatives
  oops::Parameter<int> wind_nsg = param(def.wind_nsg, this);
  // Wind inflation to compensate the Savitzky-Golay smoothing
  oops::Parameter<double> wind_inflation = param(def.wind_inflation, this);
};

// -----------------------------------------------------------------------------

// BUMP parameters
class BUMPParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BUMPParameters, oops::Parameters)

 public:
  // Internal parameters

  // General parameters
  oops::Parameter<GeneralSection> general{"general", GeneralSection(), this};
  // I/O parameters
  oops::Parameter<IOSection> io{"io", IOSection(), this};
  // Drivers parameters
  oops::Parameter<DriversSection> drivers{"drivers", DriversSection(), this};
  // Model parameters
  oops::Parameter<ModelSection> model{"model", ModelSection(), this};
  // Ensemble sizes parameters
  oops::Parameter<EnsembleSizesSection> ensembleSizes{"ensemble sizes", EnsembleSizesSection(),
    this};
  // Sampling parameters
  oops::Parameter<SamplingSection> sampling{"sampling", SamplingSection(), this};
  // Diagnostics parameters
  oops::Parameter<DiagnosticsSection> diagnostics{"diagnostics", DiagnosticsSection(), this};
  // Vertical balance parameters
  oops::Parameter<VerticalBalanceSection> verticalBalance{"vertical balance",
    VerticalBalanceSection(), this};
  // Variance parameters
  oops::Parameter<VarianceSection> variance{"variance", VarianceSection(), this};
  // Optimality test parameters
  oops::Parameter<OptimalityTestSection> optimalityTest{"optimality test", OptimalityTestSection(),
    this};
  // Fit parameters
  oops::Parameter<FitSection> fit{"fit", FitSection(), this};
  // Local profiles parameters
  oops::Parameter<std::vector<LocalProfileParameters>> localProfiles{"local profiles", {}, this};
  // NICAS parameters
  oops::Parameter<NICASSection> nicas{"nicas", NICASSection(), this};
  // Psichitouv parameters
  oops::Parameter<PsichitouvSection> psichitouv{"psichitouv", PsichitouvSection(), this};
  // Dirac parameters
  oops::Parameter<std::vector<DiracPointParameters>> dirac{"dirac", {}, this};

  // External parameters

  // Missing real value
  oops::Parameter<double> msvalr{"msvalr", util::missingValue<double>(), this};
  // Grids
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> grids{"grids", this};
  // Input ATLAS files
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> inputAtlasFilesConf{
    "input atlas files", this};
  // Input model files
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> inputModelFilesConf{
    "input model files", this};
  // Output ATLAS files
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> outputAtlasFilesConf{
    "output atlas files", this};
  // Output model files
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> outputModelFilesConf{
    "output model files", this};
};

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
