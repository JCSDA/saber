/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <math.h>

#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"

#include "oops/util/missingValues.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

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

// Specific off-diagonal weight elemental parameters
class SpecWgtParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SpecWgtParameters, oops::Parameters)

 public:
  // Variables pair
  oops::RequiredParameter<std::vector<std::string>> variablesPair{"variables pair", this};
  // Weight
  oops::RequiredParameter<double> weight{"weight", this};
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
  oops::OptionalParameter<int> level{"level", this};
  // Diracs variable indices
  oops::RequiredParameter<std::string> variable{"variable", this};
};

// -----------------------------------------------------------------------------
// BUMP parameters sections
// -----------------------------------------------------------------------------

// General section
class GeneralSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeneralSection, oops::Parameters)

 public:
  // Add colors to the log (for display on terminal)
  oops::Parameter<bool> color_log{"color log", false, this};
  // Stream test messages into a dedicated channel
  oops::Parameter<bool> testing{"testing", false, this};
  // Default seed for random numbers (0 for time-dependent seed)
  oops::Parameter<int> default_seed{"default seed", 140587, this};
  // Inter-compilers reproducibility
  oops::Parameter<bool> repro_ops{"reproducibility operators", true, this};
  // Reproducibility threshold
  oops::Parameter<double> repro_th{"reproducibility threshold", 1.0e-12, this};
  // Use timers
  oops::Parameter<bool> timers{"timers", false, this};
  // Universe radius [in meters]
  oops::Parameter<double> universe_radius{"universe length-scale", 6371229*M_PI, this};
  // Use deprecated hull (for backward compatibility)
  oops::Parameter<bool> deprecated_hull{"deprecated hull", true, this};
  // Use deprecated work grid (for backward compatibility)
  oops::Parameter<bool> deprecated_work_grid{"deprecated work grid", false, this};
};

// -----------------------------------------------------------------------------

// I/O section
class IOSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(IOSection, oops::Parameters)

 public:
  // Data directory
  oops::Parameter<std::string> data_directory{"data directory", ".", this};
  // Data prefix
  oops::Parameter<std::string> files_prefix{"files prefix", "", this};
  // Write in new files
  oops::Parameter<bool> new_files{"new files", true, this};
  // Parallel NetCDF I/O
  oops::Parameter<bool> parallel_netcdf{"parallel netcdf", true, this};
  // Number of I/O processors
  oops::Parameter<int> nprocio{"io tasks", 20, this};
  // Alias
  oops::Parameter<std::vector<AliasParameters>> alias{"alias", {}, this};
  // Sampling file
  oops::Parameter<std::string> fname_samp{"overriding sampling file", "", this};
  // Vertical covariance files
  oops::Parameter<std::vector<std::string>> fname_vbal_cov{"overriding vertical covariance file",
    {}, this};
  // Vertical balance file
  oops::Parameter<std::string> fname_vbal{"overriding vertical balance file", "", this};
  // Ensemble 1 moments files
  oops::Parameter<std::vector<std::string>> fname_mom{"overriding moments file", {}, this};
  // Averaged statistics file
  oops::Parameter<std::string> fname_avg{"overriding averaged statistics file", "", this};
  // Universe radius file
  oops::Parameter<std::string> fname_universe_radius{"overriding universe radius file", "", this};
  // NICAS file
  oops::Parameter<std::string> fname_nicas{"overriding nicas file", "", this};
  // Psichitouv transform file
  oops::Parameter<std::string> fname_wind{"overriding psichitouv file", "", this};
  // GSI data file
  oops::Parameter<std::string> fname_gsi_data{"gsi data file", "", this};
  // GSI namelist
  oops::Parameter<std::string> fname_gsi_nam{"gsi namelist", "", this};
};

// -----------------------------------------------------------------------------

// Drivers section
class DriversSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DriversSection, oops::Parameters)

 public:
  // Compute covariance
  oops::Parameter<bool> compute_cov{"compute covariance", false, this};
  // Compute correlation
  oops::Parameter<bool> compute_cor{"compute correlation", false, this};
  // Compute localization
  oops::Parameter<bool> compute_loc{"compute localization", false, this};
  // Compute hybrid weights
  oops::Parameter<bool> compute_hyb{"compute hybrid weights", false, this};
  // Hybrid term source ('randomized static' or 'lowres ensemble')
  oops::Parameter<std::string> hybrid_source{"hybrid source", "", this};
  // Multivariate strategy ('univariate', 'duplicated', 'duplicated and weighted' or 'crossed')
  oops::Parameter<std::string> strategy{"multivariate strategy", "", this};
  // New normality test
  oops::Parameter<bool> new_normality{"compute normality", false, this};
  // Read local sampling
  oops::Parameter<bool> load_samp_local{"read local sampling", false, this};
  // Read global sampling
  oops::Parameter<bool> load_samp_global{"read global sampling", false, this};
  // Write local sampling
  oops::Parameter<bool> write_samp_local{"write local sampling", false, this};
  // Write global sampling
  oops::Parameter<bool> write_samp_global{"write global sampling", false, this};
  // Write sampling grids
  oops::Parameter<bool> write_samp_grids{"write sampling grids", false, this};
  // New vertical covariance
  oops::Parameter<bool> new_vbal_cov{"compute vertical covariance", false, this};
  // Read local vertical covariance
  oops::Parameter<bool> load_vbal_cov{"read vertical covariance", false, this};
  // Write local vertical covariance
  oops::Parameter<bool> write_vbal_cov{"write vertical covariance", false, this};
  // Compute vertical balance operator
  oops::Parameter<bool> new_vbal{"compute vertical balance", false, this};
  // Read local vertical balance operator
  oops::Parameter<bool> load_vbal{"read vertical balance", false, this};
  // Write vertical balance operator
  oops::Parameter<bool> write_vbal{"write vertical balance", false, this};
  // Compute variance
  oops::Parameter<bool> new_var{"compute variance", false, this};
  // Compute moments
  oops::Parameter<bool> new_mom{"compute moments", false, this};
  // Read sampling moments
  oops::Parameter<bool> load_mom{"read moments", false, this};
  // Write sampling moments
  oops::Parameter<bool> write_mom{"write moments", false, this};
  // Write averaged statistics
  oops::Parameter<bool> write_avg{"write averaged statistics", false, this};
  // Write HDIAG diagnostics
  oops::Parameter<bool> write_hdiag{"write diagnostics", false, this};
  // Write HDIAG components detail
  oops::Parameter<bool> write_hdiag_detail{"write diagnostics detail", false, this};
  // Write HDIAG diagnostics in yaml file
  oops::Parameter<bool> write_hdiag_yaml{"write diagnostics in yaml", false, this};
  // Read universe radius
  oops::Parameter<bool> load_universe_radius{"read universe radius", false, this};
  // Write universe radius
  oops::Parameter<bool> write_universe_radius{"write universe radius", false, this};
  // Compute NICAS
  oops::Parameter<bool> new_nicas{"compute nicas", false, this};
  // Read local NICAS parameters
  oops::Parameter<bool> load_nicas_local{"read local nicas", false, this};
  // Read global NICAS parameters
  oops::Parameter<bool> load_nicas_global{"read global nicas", false, this};
  // Write local NICAS parameters
  oops::Parameter<bool> write_nicas_local{"write local nicas", false, this};
  // Write global NICAS parameters
  oops::Parameter<bool> write_nicas_global{"write global nicas", false, this};
  // Write NICAS grids
  oops::Parameter<bool> write_nicas_grids{"write nicas grids", false, this};
  // Write NICAS steps
  oops::Parameter<bool> write_nicas_steps{"write nicas steps", false, this};
  // Compute wind transform
  oops::Parameter<bool> new_wind{"compute psichitouv", false, this};
  // Read local wind transform
  oops::Parameter<bool> load_wind_local{"read local psichitouv", false, this};
  // Read global wind transform
  oops::Parameter<bool> load_wind_global{"read global psichitouv", false, this};
  // Write local wind transform
  oops::Parameter<bool> write_wind_local{"write local psichitouv", false, this};
  // Write global wind transform
  oops::Parameter<bool> write_wind_global{"write global psichitouv", false, this};
  // Test vertical balance inverse
  oops::Parameter<bool> check_vbal{"vertical balance inverse test", false, this};
  // Test adjoints
  oops::Parameter<bool> check_adjoints{"adjoints test", false, this};
  // Test NICAS normalization (number of tests)
  oops::Parameter<int> check_normalization{"normalization test", 0, this};
  // Test NICAS application on diracs
  oops::Parameter<bool> check_dirac{"internal dirac test", false, this};
  // Test NICAS randomization
  oops::Parameter<bool> check_randomization{"randomization test", false, this};
  // Test HDIAG-NICAS consistency
  oops::Parameter<bool> check_consistency{"internal consistency test", false, this};
  // Test HDIAG optimality
  oops::Parameter<bool> check_optimality{"localization optimality test", false, this};
  // Interpolate vertical balance, standard-deviation or length-scales from GSI data
  oops::Parameter<bool> from_gsi{"interpolate from gsi data", false, this};
};

// -----------------------------------------------------------------------------

// Model section
class ModelSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelSection, oops::Parameters)

 public:
  // Variables names
  oops::Parameter<std::vector<std::string>> variables{"variables", {}, this};
  // 2D variables names
  oops::Parameter<std::vector<std::string>> var2d{"2d variables", {}, this};
  // Nearest 3D level
  oops::OptionalParameter<std::string> nearest3dLevel{"nearest 3d level", this};
  // Groups of variables
  oops::OptionalParameter<std::vector<GroupParameters>> groups{"groups", this};
  // Check that sampling couples and interpolations do not cross mask boundaries
  oops::Parameter<bool> mask_check{"do not cross mask boundaries", false, this};
};

// -----------------------------------------------------------------------------

// Ensemble size section
class EnsembleSizesSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(EnsembleSizesSection, oops::Parameters)

 public:
  // Ensemble size
  oops::Parameter<int> ens_ne{"total ensemble size", 0, this};
  // Ensemble sub-ensembles number
  oops::Parameter<int> ens_nsub{"sub-ensembles", 1, this};
};

// -----------------------------------------------------------------------------

// Mask parameters
class MaskParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(MaskParameters, oops::Parameters)

 public:
  // Mask restriction type
  oops::RequiredParameter<std::string> mask_type{"type", this};
  // Mask threshold
  oops::Parameter<double> mask_th{"threshold", 0.0, this};
  // Mask threshold side ('lower' if mask_th is the lower bound, resp. 'upper')
  oops::Parameter<std::string> mask_lu{"side", "", this};
  // Mask variable
  oops::Parameter<std::string> mask_variable{"variable", "", this};
};

// -----------------------------------------------------------------------------

// Sampling section
class SamplingSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SamplingSection, oops::Parameters)

 public:
  // Computation grid size
  oops::Parameter<int> nc1{"computation grid size", 0, this};
  // Diagnostic grid size
  oops::Parameter<int> nc2{"diagnostic grid size", 0, this};
  // Number of distance classes
  oops::Parameter<int> nc3{"distance classes", 0, this};
  // Number of angular sectors
  oops::Parameter<int> nc4{"angular sectors", 1, this};
  // Class size (for sam_type='hor'), should be larger than the typical grid cell size [in meters]
  oops::Parameter<double> dc{"distance class width", 0.0, this};
  // Reduced number of levels for diagnostics
  oops::Parameter<int> nl0r{"reduced levels", 0, this};
  // Enable nl0r larger than 15 (large memory footprint)
  oops::Parameter<bool> enable_large_nl0r{"enable large number of reduced levels", false, this};
  // Activate local diagnostics
  oops::Parameter<bool> local_diag{"local diagnostic", false, this};
  // Local diagnostics calculation radius [in meters]
  oops::Parameter<double> local_rad{"averaging length-scale", 0.0, this};
  // Local diagnostics calculation latitude band half-width [in degrees]
  oops::Parameter<double> local_dlat{"averaging latitude width", 0.0, this};
  // Maximum number of random number draws
  oops::Parameter<int> irmax{"max number of draws", 10000000, this};
  // Vertical balance C2B to C0A interpolation type ('c0': C0 mesh-based, 'c1': C1 mesh-based
  // or 'si': smooth interpolation)
  oops::Parameter<std::string> interp_type{"interpolation type", "c0", this};
  // Sampling masks
  oops::Parameter<std::vector<MaskParameters>> masks{"masks", {}, this};
  // Threshold on vertically contiguous points for the mask (0 to skip the test)
  oops::Parameter<int> ncontig_th{"contiguous levels threshold", 0, this};
};

// -----------------------------------------------------------------------------

// Diagnostics section
class DiagnosticsSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DiagnosticsSection, oops::Parameters)

 public:
  // Ensemble size
  oops::Parameter<int> ne{"target ensemble size", 0, this};
  // Gaussian approximation for asymptotic quantities
  oops::Parameter<bool> gau_approx{"gaussian approximation", false, this};
  // Localization option ('default', 'from_squared_correlation', 'nice_with_table' and
  // 'nice_without_table')
  oops::Parameter<std::string> loc_option{"localization option", "default", this};
  // Threshold on generalized kurtosis (3.0 = Gaussian distribution)
  oops::Parameter<double> gen_kurt_th{"generalized kurtosis threshold",
    std::numeric_limits<double>().max(), this};
  // Number of bins for averaged statistics histograms
  oops::Parameter<int> avg_nbins{"histogram bins", 0, this};
  // Support radius scaling in CMAT from HDIAG
  oops::Parameter<double> lengths_scaling{"diagnosed lengths scaling", 1.0, this};
};

// -----------------------------------------------------------------------------

// Vertical balance block parameters
class VerticalBalanceBlockParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VerticalBalanceBlockParameters, oops::Parameters)

 public:
  // Balanced variable
  oops::RequiredParameter<std::string> balanced{"balanced variable", this};
  // Unbalanced variable
  oops::RequiredParameter<std::string> unbalanced{"unbalanced variable", this};
  // Diagonal auto-covariance for the inversion
  oops::Parameter<bool> diag_auto{"diagonal autocovariance", false, this};
  // Diagonal regression
  oops::Parameter<bool> diag_reg{"diagonal regression", false, this};
  // Scalar coefficients for identity vertical balance
  oops::Parameter<double> id_coef{"identity block weight", 1.0, this};
};

// -----------------------------------------------------------------------------

// Vertical balance section
class VerticalBalanceSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VerticalBalanceSection, oops::Parameters)

 public:
  // Vertical balance parameters
  oops::Parameter<std::vector<VerticalBalanceBlockParameters>> vbal{"vbal", {}, this};
  // Pseudo-inverse for auto-covariance
  oops::Parameter<bool> vbal_pseudo_inv{"pseudo inverse", false, this};
  // Dominant mode for pseudo-inverse
  oops::Parameter<int> vbal_pseudo_inv_mmax{"dominant mode", 0, this};
  // Variance threshold to compute the dominant mode for pseudo-inverse
  oops::Parameter<double> vbal_pseudo_inv_var_th{"variance threshold", 0.0, this};
  // Identity vertical balance for tests
  oops::Parameter<bool> vbal_id{"identity blocks", false, this};
};

// -----------------------------------------------------------------------------

// Variance section
class VarianceSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VarianceSection, oops::Parameters)

 public:
  // Force specific variance
  oops::Parameter<bool> forced_var{"explicit stddev", false, this};
  // Forced standard-deviation
  oops::Parameter<std::vector<VarsValueOrProfileParameters>> stddev{"stddev", {}, this};
  // Filter variance
  oops::Parameter<bool> var_filter{"objective filtering", false, this};
  // Number of iterations for the variance filtering (0 for uniform variance)
  oops::Parameter<int> var_niter{"filtering iterations", -1, this};
  // Number of passes for the variance filtering (0 for uniform variance)
  oops::Parameter<int> var_npass{"filtering passes", -1, this};
  // Variance initial filtering support radius [in meters]
  oops::Parameter<std::vector<VarsValueOrProfileParameters>> var_rhflt{"initial length-scale", {},
    this};
};

// -----------------------------------------------------------------------------

// Optimality test section
class OptimalityTestSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OptimalityTestSection, oops::Parameters)

 public:
  // Number of length-scale factors for optimization
  oops::Parameter<int> optimality_nfac{"half number of factors", 1, this};
  // Increments of length-scale factors for optimization
  oops::Parameter<double> optimality_delta{"factors increment", 0.05, this};
  // Number of test vectors for optimization
  oops::Parameter<int> optimality_ntest{"test vectors", 10, this};
};

// -----------------------------------------------------------------------------

// Fit section
class FitSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FitSection, oops::Parameters)

 public:
  // Threshold to filter out lower raw values
  oops::Parameter<double> diag_raw_th{"raw diagnostic lower threshold",
    -std::numeric_limits<double>().max(), this};
  // Horizontal filtering suport radius [in meters]
  oops::Parameter<double> diag_rhflt{"horizontal filtering length-scale", 0.0, this};
  // Vertical filtering support radius
  oops::Parameter<double> diag_rvflt{"vertical filtering length-scale", 0.0, this};
  // Number of levels between interpolation levels
  oops::Parameter<int> fit_dl0{"vertical stride", 1, this};
  // Number of components in the fit function
  oops::Parameter<int> fit_ncmp{"number of components", 1, this};
};

// -----------------------------------------------------------------------------

// NICAS section
class NICASSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(NICASSection, oops::Parameters)

 public:
  // Resolution
  oops::Parameter<double> resol{"resolution", 0.0, this};
  // Filter mode
  oops::Parameter<bool> filter_mode{"filter mode", false, this};
  // Resolution for the NICAS filter
  oops::Parameter<double> filter_resol{"filter resolution", 5.0, this};
  // NICAS draw type ('random' or 'regular')
  oops::Parameter<std::string> nicas_draw_type{"grid type", "regular", this};
  // Force specific support radii
  oops::Parameter<bool> forced_radii{"explicit length-scales", false, this};
  // Forced horizontal support radius [in meters]
  oops::Parameter<std::vector<GroupsValueOrProfileParameters>> rh{"horizontal length-scale", {},
    this};
  // Forced vertical support radius
  oops::Parameter<std::vector<GroupsValueOrProfileParameters>> rv{"vertical length-scale", {},
    this};
  // Force specific support radii for sampling
  oops::Parameter<bool> forced_radii_sampling{"explicit length-scales for sampling", false, this};
  // Forced horizontal support radius for sampling [in meters]
  oops::Parameter<std::vector<GroupsValueOrProfileParameters>>
    rhs{"horizontal length-scale for sampling", {}, this};
  // Forced vertical support radius for sampling
  oops::Parameter<std::vector<GroupsValueOrProfileParameters>>
    rvs{"vertical length-scale for sampling", {}, this};
  // Default off-diagonal weight
  oops::Parameter<double> defaultWeight{"default off-diagonal weight", 0.0, this};
  // Specific off-diagonal weights
  oops::Parameter<std::vector<SpecWgtParameters>> specWeights{"specific off-diagonal weights", {},
    this};
  // NICAS C1B to C0A default interpolation type ('c0': C0 mesh-based, 'c1': C1 mesh-based
  // or 'si': smooth interpolation)
  oops::Parameter<std::string> default_interp_type{"default interpolation type", "c0", this};
  // NICAS C1B to C0A interpolation type ('c0': C0 mesh-based, 'c1': C1 mesh-based
  // or 'si': smooth interpolation)
  oops::Parameter<std::vector<GroupsTypeParameters>> interp_type{"interpolation type", {},
    this};
  // Factor to get interpolation radius from convolution radius if nicas_interp_type = 'si'
  oops::Parameter<double> nicas_si_factor{"smooth interpolation factor", 0.8, this};
  // Normalization randomization size
  oops::Parameter<int> norm_rand_size{"normalization randomization size", 0, this};
  // Horizontal NICAS interpolation test
  oops::Parameter<bool> interp_test{"horizontal interpolation test", false, this};
  // Overriding component in file
  oops::Parameter<int> file_component{"overriding component in file", 0, this};
  // Same horizontal convolution for all levels, no vertical convolution
  oops::Parameter<bool> same_horizontal{"same horizontal convolution", false, this};
  // Similar levels threshold (relative tolerance on horizontal length-scale)
  oops::Parameter<double> sim_levs_th{"similar levels threshold", 0.0, this};
  // Read/write interpolation in global file
  oops::Parameter<bool> interp_in_global_file{"interpolation in global file", false, this};
  // Number of runs for setup timing
  oops::Parameter<int> nicas_setup_timings{"setup timings", 1, this};
  // Number of runs for application timing
  oops::Parameter<int> nicas_application_timings{"application timings", 1, this};
};

// -----------------------------------------------------------------------------

// Psichitouv section
class PsichitouvSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(PsichitouvSection, oops::Parameters)

 public:
  // Dipole test (bypass the adjoint)
  oops::Parameter<bool> wind_dipole_test{"dipole test", false, this};
};

// -----------------------------------------------------------------------------

// External section
class ExternalSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ExternalSection, oops::Parameters)

 public:
  // Missing real value
  oops::Parameter<double> msvalr{"msvalr", util::missingValue<double>(), this};
  // Iterative algorithm (ensemble members loaded sequentially)
  oops::Parameter<bool> iterative_algo{"iterative algorithm", false, this};
  // Vertical coordinate field name in geometry fields
  oops::Parameter<std::string> vert_coord_name{"vertical coordinate name", "", this};
  // Geographical mask name in geometry fields
  oops::Parameter<std::string> gmask_name{"geographical mask name", "", this};
};

// -----------------------------------------------------------------------------

// BUMP parameters
class BUMPParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BUMPParameters, oops::Parameters)

 public:
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
  oops::Parameter<ExternalSection> external{"external", ExternalSection(), this};

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

extern "C" {
  void bump_config_init_f90(eckit::LocalConfiguration *);
}

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
