/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_BUMP_BUMP_H_
#define SABER_BUMP_BUMP_H_

#include <omp.h>

#include <algorithm>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/config/Configuration.h"

#include "oops/base/Increment.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

#include "saber/bump/type_bump.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------
/// Parameters describing BUMP components input (from generic text files).
class BUMPInputNcmpParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BUMPInputNcmpParameters, oops::Parameters)

 public:
  /// File path
  oops::RequiredParameter<std::string> filepath{"filepath", this};
};

// -----------------------------------------------------------------------------
/// Parameters describing BUMP parameters input (from model Increment files).
template <typename MODEL> class BUMPInputParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BUMPInputParameters, oops::Parameters)
  typedef typename oops::Increment<MODEL>::ReadParameters_ ReadParameters_;

 public:
  /// Parameter name.
  oops::RequiredParameter<std::string> param{"parameter", this};
  /// Component index
  oops::Parameter<int> component{"component", 1, this};
  /// Parameters used for reading Increment.
  ReadParameters_ incread{this};
};

// -----------------------------------------------------------------------------
/// Parameters describing BUMP components output (to generic text files)
class BUMPOutputNcmpParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BUMPOutputNcmpParameters, oops::Parameters)

 public:
  /// File path
  oops::RequiredParameter<std::string> filepath{"filepath", this};
};

// -----------------------------------------------------------------------------
/// Parameters describing BUMP parameters output (to model Increment files)
template <typename MODEL> class BUMPOutputParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BUMPOutputParameters, oops::Parameters)
  typedef typename oops::Increment<MODEL>::WriteParameters_ WriteParameters_;

 public:
  /// Parameter name.
  oops::RequiredParameter<std::string> param{"parameter", this};
  /// Component index
  oops::Parameter<int> component{"component", 1, this};
  /// Parameters used for writing Increment.
  WriteParameters_ incwrite{this};
};

template <typename MODEL> class BUMP_Parameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BUMP_Parameters, oops::Parameters)

 public:
  // Internal parameters

  // general_param

  // Data directory
  oops::OptionalParameter<std::string> datadir{"datadir", this};
  // Files prefix
  oops::OptionalParameter<std::string> prefix{"prefix", this};
  // Model name ('aro', 'arp', 'fv3', 'gem', 'geos', 'gfs', 'ifs', 'mpas', 'nemo', 'norcpm',
  // 'online', 'qg, 'res' or 'wrf')
  oops::OptionalParameter<std::string> model{"model", this};
  // Verbosity level ('all', 'main' or 'none')
  oops::OptionalParameter<std::string> verbosity{"verbosity", this};
  // Add colors to the log (for display on terminal)
  oops::OptionalParameter<bool> colorlog{"colorlog", this};
  // Stream test messages into a dedicated channel
  oops::OptionalParameter<bool> testing{"testing", this};
  // Default seed for random numbers
  oops::OptionalParameter<bool> default_seed{"default_seed", this};
  // Inter-compilers reproducibility
  oops::OptionalParameter<bool> repro{"repro", this};
  // Reproducibility threshold
  oops::OptionalParameter<double> rth{"rth", this};
  // Parallel NetCDF I/O
  oops::OptionalParameter<bool> parallel_io{"parallel_io", this};
  // Number of I/O processors
  oops::OptionalParameter<int> nprocio{"nprocio", this};
  // Universe radius [in meters]
  oops::OptionalParameter<double> universe_rad{"universe_rad", this};
  // Write subset Sc0 fields (full grid) using BUMP I/O
  oops::OptionalParameter<bool> write_c0{"write_c0", this};

  // driver_param

  // Localization/hybridization to compute ('cor', 'loc', 'hyb-rnd' or 'hyb-ens')
  oops::OptionalParameter<std::string> method{"method", this};
  // Localization strategy ('diag_all', 'common', 'common_weighted', 'specific_univariate' or
  // 'specific_multivariate')
  oops::OptionalParameter<std::string> strategy{"strategy", this};
  // New normality test
  oops::OptionalParameter<bool> new_normality{"new_normality", this};
  // New vertical covariance
  oops::OptionalParameter<bool> new_vbal_cov{"new_vbal_cov", this};
  // Update vertical covariance sequentially
  oops::OptionalParameter<bool> update_vbal_cov{"update_vbal_cov", this};
  // Load local vertical covariance
  oops::OptionalParameter<bool> load_vbal_cov{"load_vbal_cov", this};
  // Write local vertical covariancee
  oops::OptionalParameter<bool> write_vbal_cov{"write_vbal_cov", this};
  // Compute new vertical balance operator
  oops::OptionalParameter<bool> new_vbal{"new_vbal", this};
  // Load local vertical balance operator
  oops::OptionalParameter<bool> load_vbal{"load_vbal", this};
  // Write vertical balance operator
  oops::OptionalParameter<bool> write_vbal{"write_vbal", this};
  // Compute new variance
  oops::OptionalParameter<bool> new_var{"new_var", this};
  // Update variance sequentially
  oops::OptionalParameter<bool> update_var{"update_var", this};
  // Load variance
  oops::OptionalParameter<bool> load_var{"load_var", this};
  // Write variance
  oops::OptionalParameter<bool> write_var{"write_var", this};
  // Compute new sampling moments
  oops::OptionalParameter<bool> new_mom{"new_mom", this};
  // Update sampling moments sequentially
  oops::OptionalParameter<bool> update_mom{"update_mom", this};
  // Load sampling moments
  oops::OptionalParameter<bool> load_mom{"load_mom", this};
  // Write sampling moments
  oops::OptionalParameter<bool> write_mom{"write_mom", this};
  // Compute new HDIAG diagnostics
  oops::OptionalParameter<bool> new_hdiag{"new_hdiag", this};
  // Write HDIAG diagnostics
  oops::OptionalParameter<bool> write_hdiag{"write_hdiag", this};
  // Compute new NICAS parameters
  oops::OptionalParameter<bool> new_nicas{"new_nicas", this};
  // Load local NICAS parameters
  oops::OptionalParameter<bool> load_nicas_local{"load_nicas_local", this};
  // Load global NICAS parameters
  oops::OptionalParameter<bool> load_nicas_global{"load_nicas_global", this};
  // Write local NICAS parameters
  oops::OptionalParameter<bool> write_nicas_local{"write_nicas_local", this};
  // Write global NICAS parameters
  oops::OptionalParameter<bool> write_nicas_global{"write_nicas_global", this};
  // Compute wind transform
  oops::OptionalParameter<bool> new_wind{"new_wind", this};
  // Load local wind transform
  oops::OptionalParameter<bool> load_wind_local{"load_wind_local", this};
  // Write local wind transform
  oops::OptionalParameter<bool> write_wind_local{"write_wind_local", this};
  // Test vertical balance inverse and adjoint
  oops::OptionalParameter<bool> check_vbal{"check_vbal", this};
  // Test NICAS adjoints
  oops::OptionalParameter<bool> check_adjoints{"check_adjoints", this};
  // Test NICAS normalization (number of tests)
  oops::OptionalParameter<int> check_normalization{"check_normalization", this};
  // Test NICAS application on diracs
  oops::OptionalParameter<bool> check_dirac{"check_dirac", this};
  // Test NICAS randomization
  oops::OptionalParameter<bool> check_randomization{"check_randomization", this};
  // Test HDIAG-NICAS consistency
  oops::OptionalParameter<bool> check_consistency{"check_consistency", this};
  // Test HDIAG optimality
  oops::OptionalParameter<bool> check_optimality{"check_optimality", this};
  // Test BUMP with no grid point on the last MPI task
  oops::OptionalParameter<bool> check_no_point_mpi{"check_no_point_mpi", this};
  // Test BUMP with all grid points masked on half of the domain
  oops::OptionalParameter<bool> check_no_point_mask{"check_no_point_mask", this};
  // Test set_parameter interface
  oops::OptionalParameter<bool> check_set_param{"check_set_param", this};
  // Test get_parameter interface
  oops::OptionalParameter<bool> check_get_param{"check_get_param", this};
  // Test apply_vbal interfaces
  oops::OptionalParameter<bool> check_apply_vbal{"check_apply_vbal", this};
  // Test apply_stddev interfaces
  oops::OptionalParameter<bool> check_apply_stddev{"check_apply_stddev", this};
  // Test apply_nicas interfaces
  oops::OptionalParameter<bool> check_apply_nicas{"check_apply_nicas", this};

  // files_param

  // Variance files
  oops::OptionalParameter<std::vector<std::string>> fname_var{"fname_var", this};
  // Sampling file
  oops::OptionalParameter<std::string> fname_samp{"fname_samp", this};
  // Vertical covariance files
  oops::OptionalParameter<std::vector<std::string>> fname_vbal_cov{"fname_vbal_cov", this};
  // Vertical balance file
  oops::OptionalParameter<std::string> fname_vbal{"fname_vbal", this};
  // Moments files
  oops::OptionalParameter<std::vector<std::string>> fname_mom{"fname_mom", this};
  // NICAS file
  oops::OptionalParameter<std::string> fname_nicas{"fname_nicas", this};
  // Wind transform file
  oops::OptionalParameter<std::string> fname_wind{"fname_wind", this};

  // model_param

  // Number of levels
  oops::OptionalParameter<int> nl0{"nl0", this};
  // Levels
  oops::OptionalParameter<std::vector<int>> levs{"levs", this};
  // Level for 2D variables ('first' or 'last')
  oops::OptionalParameter<std::string> lev2d{"lev2d", this};
  // Use pressure logarithm as vertical coordinate (model level if .false.)
  oops::OptionalParameter<bool> logpres{"logpres", this};
  // Number of variables
  oops::OptionalParameter<int> nv{"nv", this};
  // Variables names
  oops::OptionalParameter<std::vector<std::string>> variables{"variables", this};
  // Variable change
  oops::OptionalParameter<std::string> variable_change{"variable_change", this};
  // Do not use geometry mask
  oops::OptionalParameter<bool> nomask{"nomask", this};
  // I/O keys
  oops::OptionalParameter<std::vector<std::string>> io_keys{"io_keys", this};
  // I/O values
  oops::OptionalParameter<std::vector<std::string>> io_values{"io_values", this};
  // Regional domain configuration for the QG model
  oops::OptionalParameter<bool> qg_regional{"qg_regional", this};
  // Urban domain configuration for the QG model
  oops::OptionalParameter<bool> qg_urban{"qg_urban", this};

  // ens1_param

  // Ensemble 1 size
  oops::OptionalParameter<int> ens1_ne{"ens1_ne", this};
  // Ensemble 1 sub-ensembles number
  oops::OptionalParameter<int> ens1_nsub{"ens1_nsub", this};

  // ens2_param

  // Ensemble 2 size
  oops::OptionalParameter<int> ens2_ne{"ens2_ne", this};
  // Ensemble 2 sub-ensembles number
  oops::OptionalParameter<int> ens2_nsub{"ens2_nsub", this};

  // sampling_param

  // Load local sampling
  oops::OptionalParameter<bool> load_samp_local{"load_samp_local", this};
  // Load global sampling
  oops::OptionalParameter<bool> load_samp_global{"load_samp_global", this};
  // Write local sampling
  oops::OptionalParameter<bool> write_samp_local{"write_samp_local", this};
  // Write global sampling
  oops::OptionalParameter<bool> write_samp_global{"write_samp_global", this};
  // Write sampling grids
  oops::OptionalParameter<bool> write_samp_grids{"write_samp_grids", this};
  // Mask restriction type
  oops::OptionalParameter<std::string> mask_type{"mask_type", this};
  // Mask threshold side ('lower' if mask_th is the lower bound, 'upper' if mask_th is the
  // upper bound)
  oops::OptionalParameter<std::vector<std::string>> mask_lu{"mask_lu", this};
  // Mask threshold
  oops::OptionalParameter<std::vector<double>> mask_th{"mask_th", this};
  // Threshold on vertically contiguous points for sampling mask (0 to skip the test)
  oops::OptionalParameter<int> ncontig_th{"ncontig_th", this};
  // Check that sampling couples and interpolations do not cross mask boundaries
  oops::OptionalParameter<bool> mask_check{"mask_check", this};
  // Sampling draw type ('random_uniform','random_coast' or 'octahedral')
  oops::OptionalParameter<std::string> diag_draw_type{"diag_draw_type", this};
  // Length-scale to increase sampling density along coasts [in meters]
  oops::OptionalParameter<double> Lcoast{"Lcoast", this};
  // Minimum value to increase sampling density along coasts
  oops::OptionalParameter<double> rcoast{"rcoast", this};
  // Number of sampling points
  oops::OptionalParameter<int> nc1{"nc1", this};
  // Number of diagnostic points
  oops::OptionalParameter<int> nc2{"nc2", this};
  // Number of horizontal classes
  oops::OptionalParameter<int> nc3{"nc3", this};
  // Number of angular sectors
  oops::OptionalParameter<int> nc4{"nc4", this};
  // Class size (for sam_type='hor'), should be larger than the typical grid cell size [in meters]
  oops::OptionalParameter<double> dc{"dc", this};
  // Reduced number of levels for diagnostics
  oops::OptionalParameter<int> nl0r{"nl0r", this};
  // Maximum number of random number draws
  oops::OptionalParameter<int> irmax{"irmax", this};

  // diag_param

  // Ensemble size
  oops::OptionalParameter<int> ne{"ne", this};
  // Ensemble size of the hybrid term
  oops::OptionalParameter<int> ne_lr{"ne_lr", this};
  // Threshold on generalized kurtosis (3.0 = Gaussian distribution)
  oops::OptionalParameter<double> gen_kurt_th{"gen_kurt_th", this};
  // Gaussian approximation for asymptotic quantities
  oops::OptionalParameter<bool> gau_approx{"gau_approx", this};
  // Number of bins for averaged statistics histograms
  oops::OptionalParameter<int> avg_nbins{"avg_nbins", this};
  // Activation of vertical balance (ordered line by line in the lower triangular formulation)
  oops::OptionalParameter<std::vector<bool>> vbal_block{"vbal_block", this};
  // Vertical balance diagnostic radius [in meters]
  oops::OptionalParameter<double> vbal_rad{"vbal_rad", this};
  // Vertical balance diagnostic latitude band half-width [in degrees]
  oops::OptionalParameter<double> vbal_dlat{"vbal_dlat", this};
  // Diagonal auto-covariance for the inversion
  oops::OptionalParameter<std::vector<bool>> vbal_diag_auto{"vbal_diag_auto", this};
  // Diagonal regression
  oops::OptionalParameter<std::vector<bool>> vbal_diag_reg{"vbal_diag_reg", this};
  // Pseudo-inverse for auto-covariance
  oops::OptionalParameter<bool> vbal_pseudo_inv{"vbal_pseudo_inv", this};
  // Dominant mode for pseudo-inverse
  oops::OptionalParameter<int> vbal_pseudo_inv_mmax{"vbal_pseudo_inv_mmax", this};
  // Variance threshold to compute the dominant mode for pseudo-inverse
  oops::OptionalParameter<double> vbal_pseudo_inv_var_th{"vbal_pseudo_inv_var_th", this};
  // Force specific variance
  oops::OptionalParameter<bool> forced_var{"forced_var", this};
  // Forced standard-deviation
  oops::OptionalParameter<eckit::LocalConfiguration> stddev{"stddev", this};
  // Filter variance
  oops::OptionalParameter<bool> var_filter{"var_filter", this};
  // Number of iterations for the variance filtering (0 for uniform variance)
  oops::OptionalParameter<int> var_niter{"var_niter", this};
  // Number of passes for the variance filtering (0 for uniform variance)
  oops::OptionalParameter<int> var_npass{"var_npass", this};
  // Variance initial filtering support radius [in meters]
  oops::OptionalParameter<eckit::LocalConfiguration> var_rhflt{"var_rhflt", this};
  // Activate local diagnostics
  oops::OptionalParameter<bool> local_diag{"local_diag", this};
  // Local diagnostics calculation radius [in meters]
  oops::OptionalParameter<double> local_rad{"local_rad", this};
  // Local diagnostics calculation latitude band half-width [in degrees]
  oops::OptionalParameter<double> local_dlat{"local_dlat", this};

  // fit_param

  // Horizontal filtering suport radius [in meters]
  oops::OptionalParameter<double> diag_rhflt{"diag_rhflt", this};
  // Vertical filtering support radius
  oops::OptionalParameter<double> diag_rvflt{"diag_rvflt", this};
  // Number of levels between interpolation levels
  oops::OptionalParameter<int> fit_dl0{"fit_dl0", this};
  // Number of components in the fit function
  oops::OptionalParameter<eckit::LocalConfiguration> fit_ncmp{"fit_ncmp", this};
  // Write HDIAG components detail
  oops::OptionalParameter<bool> write_hdiag_detail{"write_hdiag_detail", this};
  // Update localization for hybridization
  oops::OptionalParameter<bool> hybrid_loc_update{"hybrid_loc_update", this};

  // nicas_param

  // Resolution
  oops::OptionalParameter<double> resol{"resol", this};
  // Maximum size of the Sc1 subset
  oops::OptionalParameter<int> nc1max{"nc1max", this};
  // Subsampling draw type ('random_uniform','random_coast' or 'octahedral')
  oops::OptionalParameter<std::string> nicas_draw_type{"nicas_draw_type", this};
  // Network-base convolution calculation (distance-based if false)
  oops::OptionalParameter<bool> network{"network", this};
  // Force specific support radii
  oops::OptionalParameter<bool> forced_radii{"forced_radii", this};
  // Forced horizontal support radius [in meters]
  oops::OptionalParameter<eckit::LocalConfiguration> rh{"rh", this};
  // Forced vertical support radius
  oops::OptionalParameter<eckit::LocalConfiguration> rv{"rv", this};
  // Forced localization weights
  oops::OptionalParameter<eckit::LocalConfiguration> loc_wgt{"loc_wgt", this};
  // Minimum level
  oops::OptionalParameter<eckit::LocalConfiguration> min_lev{"min_lev", this};
  // Maximum level
  oops::OptionalParameter<eckit::LocalConfiguration> max_lev{"max_lev", this};
  // Positive-definiteness test
  oops::OptionalParameter<bool> pos_def_test{"pos_def_test", this};
  // Write NICAS grids
  oops::OptionalParameter<bool> write_nicas_grids{"write_nicas_grids", this};

  // dirac_param

  // Number of Diracs
  oops::OptionalParameter<int> ndir{"ndir", this};
  // Diracs longitudes [in degrees]
  oops::OptionalParameter<std::vector<double>> londir{"londir", this};
  // Diracs latitudes [in degrees]
  oops::OptionalParameter<std::vector<double>> latdir{"latdir", this};
  // Diracs level
  oops::OptionalParameter<std::vector<int>> levdir{"levdir", this};
  // Diracs variable indices
  oops::OptionalParameter<std::vector<int>> ivdir{"ivdir", this};

  // output_param

  // Number of neighbors for the full grid smoother
  oops::OptionalParameter<int> full_grid_smoother_nn{"full_grid_smoother_nn", this};
  // Number of local diagnostics profiles to write (for local_diag = .true.)
  oops::OptionalParameter<int> nldwv{"nldwv", this};
  // Index on model grid of the local diagnostics profiles to write
  oops::OptionalParameter<std::vector<int>> img_ldwv{"img_ldwv", this};
  // Longitudes of the local diagnostics profiles to write [in degrees]
  oops::OptionalParameter<std::vector<double>> lon_ldwv{"lon_ldwv", this};
  // Latitudes of the local diagnostics profiles to write [in degrees]
  oops::OptionalParameter<std::vector<double>> lat_ldwv{"lat_ldwv", this};
  // Name of the local diagnostics profiles to write
  oops::OptionalParameter<std::vector<std::string>> name_ldwv{"name_ldwv", this};

  // wind_param

  // Streamfunction variable name
  oops::OptionalParameter<std::string> wind_streamfunction{"wind_streamfunction", this};
  // Velocity potential variable name
  oops::OptionalParameter<std::string> wind_velocity_potential{"wind_velocity_potential", this};
  // Zonal wind variable name
  oops::OptionalParameter<std::string> wind_zonal{"wind_zonal", this};
  // Meridional variable name
  oops::OptionalParameter<std::string> wind_meridional{"wind_meridional", this};
  // Number of longitudes for the regular grid
  oops::OptionalParameter<int> wind_nlon{"wind_nlon", this};
  // Number of latitudes for the regular grid
  oops::OptionalParameter<int> wind_nlat{"wind_nlat", this};
  // Half-width of the Savitzky-Golay to compute derivatives
  oops::OptionalParameter<int> wind_nsg{"wind_nsg", this};
  // Wind inflation to compensate the Savitzky-Golay smoothing
  oops::OptionalParameter<double> wind_inflation{"wind_inflation", this};

  // External parameters

  // Universe radius (increment)
  oops::OptionalParameter<eckit::LocalConfiguration> universeRadius{"universe radius", this};
  // Components input parameters
  oops::OptionalParameter<BUMPInputNcmpParameters> inputNcmp{"input number of components", this};
  // Input parameters
  oops::OptionalParameter<std::vector<BUMPInputParameters<MODEL>>> input{"input", this};
  // Ensemble 1 parameters
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble1{"ensemble", this};
  // Ensemble 2 parameters
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble2{"lowres ensemble", this};
  // Missing value (real)
  oops::OptionalParameter<double> msvalr{"msvalr", this};
  // Grids
  oops::OptionalParameter<eckit::LocalConfiguration> grids{"grids", this};
  // Output number of components
  oops::OptionalParameter<BUMPOutputNcmpParameters> outputNcmp{"output number of components", this};
  // Output parameters
  oops::OptionalParameter<std::vector<BUMPOutputParameters<MODEL>>> output{"output", this};
  // Operators application
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> appConfs{"operators application",
    this};
};

// -----------------------------------------------------------------------------

template<typename MODEL> class BUMP {
  typedef oops::Geometry<MODEL>                           Geometry_;
  typedef oops::Increment<MODEL>                          Increment_;
  typedef BUMP_Parameters<MODEL>                          BUMP_Parameters_;
  typedef oops::State<MODEL>                              State_;
  typedef std::shared_ptr<oops::IncrementEnsemble<MODEL>> EnsemblePtr_;

 public:
  // Constructors
  BUMP(const Geometry_ &,
       const Geometry_ &,
       const oops::Variables &,
       const BUMP_Parameters_ &,
       const State_ &,
       const State_ &,
       const EnsemblePtr_ ens1 = NULL,
       const EnsemblePtr_ ens2 = NULL);

  // Copy
  explicit BUMP(BUMP &);

  // Destructor
  ~BUMP();

  // Fortran interfaces
  void addMember(const atlas::FieldSet &, const int &, const int &) const;
  void updateVbalCov(const atlas::FieldSet &, const int &) const;
  void updateVar(const atlas::FieldSet &, const int &) const;
  void updateMom(const atlas::FieldSet &, const int &, const int &) const;
  void runDrivers() const;
  void multiplyVbal(atlas::FieldSet &) const;
  void inverseMultiplyVbal(atlas::FieldSet &) const;
  void multiplyVbalAd(atlas::FieldSet &) const;
  void inverseMultiplyVbalAd(atlas::FieldSet &) const;
  void multiplyStdDev(atlas::FieldSet &) const;
  void inverseMultiplyStdDev(atlas::FieldSet &) const;
  void randomizeNicas(atlas::FieldSet &) const;
  void multiplyNicas(atlas::FieldSet &) const;
  void multiplyPsiChiToUV(atlas::FieldSet &) const;
  void multiplyPsiChiToUVAd(atlas::FieldSet &) const;
  void getNcmp(const int &, const int &, int &) const;
  void getParameter(const std::string &, const int &, const int &, atlas::FieldSet &) const;
  void setNcmp(const int &, const int &, const int &) const;
  void setParameter(const std::string &, const int &, const atlas::FieldSet &) const;
  void partialDealloc() const;

 private:
  const oops::Variables activeVars_;
  BUMP_Parameters_ params_;
  std::vector<int> keyBUMP_;
  std::vector<oops::Variables> activeVarsPerGrid_;
  bool verbosity_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
BUMP<MODEL>::BUMP(const Geometry_ & geom1,
                  const Geometry_ & geom2,
                  const oops::Variables & activeVars,
                  const BUMP_Parameters_ & params,
                  const State_ & xb,
                  const State_ & fg,
                  const EnsemblePtr_ ens1,
                  const EnsemblePtr_ ens2)
  : activeVars_(activeVars), params_(params), keyBUMP_(), activeVarsPerGrid_() {
  oops::Log::trace() << "BUMP<MODEL>::BUMP construction starting" << std::endl;

  // BUMP verbosity level ('all', 'main' or 'none')
  const boost::optional<std::string> &verbosity = params_.verbosity.value();
  if (verbosity == boost::none) {
    verbosity_ = true;
  } else {
    if (*verbosity  == "all") {
      verbosity_ = true;
    } else if (*verbosity  == "main") {
      verbosity_ = (geom1.getComm().rank() == 0);
    } else if (*verbosity  == "none") {
      verbosity_ = false;
    } else {
      ABORT("BUMP::BUMP: wrong verbosity string");
    }
  }

  // If testing is activated, replace _MPI_ and _OMP_ patterns
  const boost::optional<bool> &testing = params_.testing.value();
  if (testing != boost::none) {
    if (testing) {
      // Convert to eckit configuration
      eckit::LocalConfiguration fullConfig;
      params_.serialize(fullConfig);

      // Get number of MPI tasks and OpenMP threads
      std::string mpi(std::to_string(geom1.getComm().size()));
      std::string omp("1");
      # pragma omp parallel
      {
          omp = std::to_string(omp_get_num_threads());
      }
      if (verbosity_) oops::Log::info() << "Info     : MPI tasks:      " << mpi << std::endl;
      if (verbosity_) oops::Log::info() << "Info     : OpenMP threads: " << omp << std::endl;

      // Replace patterns
      util::seekAndReplace(fullConfig, "_MPI_", mpi);
      util::seekAndReplace(fullConfig, "_OMP_", omp);

      // Convert back to parameters
      params_.deserialize(fullConfig);
    }
  }


  // Define base increments
  Increment_ dx1(geom1, activeVars_, xb.validTime());
  Increment_ dx2(geom2, activeVars_, xb.validTime());

  // Get ensemble 1 size if ensemble 1 is available
  int ens1_ne = 0;
  if (ens1) ens1_ne = ens1->size();
  const boost::optional<eckit::LocalConfiguration> &ensembleConfig1 = params_.ensemble1.value();
  std::vector<eckit::LocalConfiguration> membersConfig1;
  if (ensembleConfig1 != boost::none) {
    // Abort if both "members" and "members from template" are specified
    if (ensembleConfig1->has("members") && ensembleConfig1->has("members from template"))
      ABORT("BUMP::BUMP: both members and members from template are specified");

    if (ensembleConfig1->has("members")) {
      // Explicit members
      ensembleConfig1->get("members", membersConfig1);
      ens1_ne = membersConfig1.size();
    } else if (ensembleConfig1->has("members from template")) {
      // Templated members
      eckit::LocalConfiguration templateConfig;
      ensembleConfig1->get("members from template", templateConfig);
      eckit::LocalConfiguration membersTemplate;
      templateConfig.get("template", membersTemplate);
      std::string pattern;
      templateConfig.get("pattern", pattern);
      templateConfig.get("nmembers", ens1_ne);
      int start = 1;
      if (templateConfig.has("start")) {
        templateConfig.get("start", start);
      }
      std::vector<int> except;
      if (templateConfig.has("except")) {
        templateConfig.get("except", except);
      }
      int zpad = 0;
      if (templateConfig.has("zero padding")) {
        templateConfig.get("zero padding", zpad);
      }
      int count = start;
      for (int ie=0; ie < ens1_ne; ++ie) {
        while (std::count(except.begin(), except.end(), count)) {
          count += 1;
        }
        eckit::LocalConfiguration memberConfig(membersTemplate);
        util::seekAndReplace(memberConfig, pattern, count, zpad);
        membersConfig1.push_back(memberConfig);
        count += 1;
      }
    } else {
      ABORT("BUMP::BUMP: ensemble 1 not specified");
    }
  }

  // Get ensemble 2 size if ensemble 2 is available
  int ens2_ne = 0;
  if (ens2) ens2_ne = ens2->size();
  const boost::optional<eckit::LocalConfiguration> &ensembleConfig2 = params_.ensemble2.value();
  std::vector<eckit::LocalConfiguration> membersConfig2;
  if (ensembleConfig2 != boost::none) {
    // Abort if both "members" and "members from template" are specified
    if (ensembleConfig2->has("members") && ensembleConfig2->has("members from template"))
      ABORT("BUMP::BUMP: both members and members from template are specified");

    if (ensembleConfig2->has("members")) {
      // Explicit members
      ensembleConfig2->get("members", membersConfig2);
      ens2_ne = membersConfig2.size();
    } else if (ensembleConfig2->has("members from template")) {
      // Templated members
      eckit::LocalConfiguration templateConfig;
      ensembleConfig2->get("members from template", templateConfig);
      eckit::LocalConfiguration membersTemplate;
      templateConfig.get("template", membersTemplate);
      std::string pattern;
      templateConfig.get("pattern", pattern);
      templateConfig.get("nmembers", ens2_ne);
      int start = 1;
      if (templateConfig.has("start")) {
        templateConfig.get("start", start);
      }
      std::vector<int> except;
      if (templateConfig.has("except")) {
        templateConfig.get("except", except);
      }
      int zpad = 0;
      if (templateConfig.has("zero padding")) {
        templateConfig.get("zero padding", zpad);
      }
      int count = start;
      for (int ie=0; ie < ens2_ne; ++ie) {
        while (std::count(except.begin(), except.end(), count)) {
          count += 1;
        }
        eckit::LocalConfiguration memberConfig(membersTemplate);
        util::seekAndReplace(memberConfig, pattern, count, zpad);
        membersConfig2.push_back(memberConfig);
        count += 1;
      }
    } else {
      ABORT("BUMP::BUMP: ensemble 2 not specified");
    }
  }

  // Read universe size
  if (verbosity_) oops::Log::info() << "Info     : Read universe radius" << std::endl;
  atlas::FieldSet universe_rad = atlas::FieldSet();
  const boost::optional<eckit::LocalConfiguration> &universeRadius = params_.universeRadius.value();
  if (universeRadius != boost::none) {
    // Read universe radius
    dx1.read(*universeRadius);

    // Get ATLAS fieldset
    for (const auto & field : dx1.fieldSet()) {
      universe_rad.add(field);
    }
  }

  // Add ensemble sizes
  eckit::LocalConfiguration conf(params_.toConfiguration());
  if (!conf.has("ens1_ne")) conf.set("ens1_ne", ens1_ne);
  if (!conf.has("ens2_ne")) conf.set("ens2_ne", ens2_ne);

  // Add missing value
  const double msvalr = util::missingValue(double());
  conf.set("msvalr", msvalr);

  // Grids
  std::vector<eckit::LocalConfiguration> grids;

  // Get global prefix
  std::string prefix;
  if (conf.has("prefix")) {
    conf.get("prefix", prefix);
  } else {
    prefix = "bump";
  }

  // Get the grids configuration from input configuration and complete it
  if (conf.has("grids")) {
    // Get grids from input configuration
    conf.get("grids", grids);
    ASSERT(grids.size() > 0);
  } else {
    // Create one empty configuration
    eckit::LocalConfiguration emptyConf;
    grids.push_back(emptyConf);
  }

  // Check grids number
  ASSERT(grids.size() > 0);

  // Loop over grids
  for (unsigned int jgrid = 0; jgrid < grids.size(); ++jgrid) {
    // Add prefix
    if (!grids[jgrid].has("prefix")) {
      std::ostringstream ss;
      ss << std::setw(2) << std::setfill('0') << jgrid;
      grids[jgrid].set("prefix", prefix + "_" + ss.str());
    }

    // Get ATLAS variable names
    std::vector<std::string> vars_atlas;
    for (const auto & field : dx1.fieldSet()) {
      vars_atlas.push_back(field.name());
    }

    // Add input variables to the grid configuration
    std::vector<std::string> vars_str;
    if (grids[jgrid].has("variables")) {
      grids[jgrid].get("variables", vars_str);
    } else {
      vars_str = vars_atlas;
      grids[jgrid].set("variables", vars_str);
    }
    grids[jgrid].set("nv", vars_str.size());

    // Save variables for each grid
    const oops::Variables gridVars(vars_str);
    activeVarsPerGrid_.push_back(gridVars);

    // Get the required number of levels add it to the grid configuration
    int nl0 = 0;
    for (const auto & field : dx1.fieldSet()) {
      if (gridVars.has(field.name())) {
        nl0 = std::max(nl0, std::max(field.levels(), 1));
      }
    }
    grids[jgrid].set("nl0", nl0);

    // Add level index for 2D fields (first or last, first by default)
    if (!grids[jgrid].has("lev2d")) {
      grids[jgrid].set("lev2d", "first");
    }

    // Print configuration for this grid
    if (verbosity_) oops::Log::info() << "Info     : Grid " << jgrid << ": " << grids[jgrid]
      << std::endl;

    // Create BUMP instance
    if (verbosity_) oops::Log::info() << "Info     : Create BUMP instance " << jgrid << std::endl;
    int keyBUMP = 0;
    bump_create_f90(keyBUMP, &geom1.getComm(),
                    geom1.functionSpace().get(),
                    geom1.extraFields().get(),
                    conf, grids[jgrid], universe_rad.get());
    keyBUMP_.push_back(keyBUMP);

    // Second geometry
    if (ens2 || (ensembleConfig2 != boost::none)) {
      bump_second_geometry_f90(keyBUMP,
                               geom2.functionSpace().get(),
                               geom2.extraFields().get());
    }
  }

  // Add members of ensemble 1
  if (ens1) {
    if (verbosity_) oops::Log::info() << "Info     : --- Add members of ensemble 1" << std::endl;
    for (int ie = 0; ie < ens1_ne; ++ie) {
      if (verbosity_) oops::Log::info() << "Info     :       Member " << ie+1 << " / " << ens1_ne
        << std::endl;
      this->addMember((*ens1)[ie].fieldSet(), ie, 1);
    }
  }

  // Add members of ensemble 2
  if (ens2) {
    if (verbosity_) oops::Log::info() << "Info     : --- Add members of ensemble 2" << std::endl;
    for (int ie = 0; ie < ens2_ne; ++ie) {
      if (verbosity_) oops::Log::info() << "Info     :       Member " << ie+1 << " / " << ens2_ne
        << std::endl;
      this->addMember((*ens2)[ie].fieldSet(), ie, 2);
    }
  }

  // Reset parameters
  params_.validateAndDeserialize(conf);

  // Read number of components
  if (verbosity_) oops::Log::info() << "Info     :     Read number of components" << std::endl;
  const boost::optional<BUMPInputNcmpParameters> &inputNcmp = params_.inputNcmp.value();
  if (inputNcmp != boost::none) {
    // Open file
    std::ifstream infile;
    infile.open(inputNcmp->filepath.value().c_str());

    if (infile.is_open()) {
      // Read file
      std::string line;
      while (std::getline(infile, line)) {
        // Split string
        std::istringstream iss(line);
        std::vector<std::string> split(std::istream_iterator<std::string>{iss},
                                       std::istream_iterator<std::string>());
        const std::string variable(split[0]);
        const int ncmp = std::stoi(split[1]);

        // Get grid and variable index
        int igrid = -1;
        int ivar = -1;
        for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
          for (size_t jvar=0; jvar < activeVarsPerGrid_[jgrid].size(); ++jvar) {
            if (activeVarsPerGrid_[jgrid][jvar] == variable) {
              igrid = jgrid;
              ivar = jvar;
              break;
            }
          }
        }
        if (igrid == -1 || ivar == -1) {
           ABORT("BUMP::BUMP: cannot find indices for variable " + variable);
        }

        // Set parameter
        this->setNcmp(igrid, ivar, ncmp);
        if (verbosity_) oops::Log::test() << "Number of input BUMP components for " << variable
          << ": " << ncmp << std::endl;
      }

      // Close file
      infile.close();
    } else {
      ABORT("BUMP::BUMP: cannot open file");
    }
  }

  // Read parameters from files
  if (verbosity_) oops::Log::info() << "Info     :     Read parameters from files" << std::endl;
  const boost::optional<std::vector<BUMPInputParameters<MODEL>>> &input = params_.input.value();
  if (input != boost::none) {
    // Set input parameters
    for (const auto & inputParam : *input) {
      // Read increment
      dx1.read(inputParam.incread);

      // Set parameter to BUMP
      const std::string & param = inputParam.param;
      const int & component = inputParam.component;
      this->setParameter(param, component, dx1.fieldSet());
      if (verbosity_) oops::Log::test() << "Norm of input BUMP parameter " << param << " - "
        << component << ": "<< dx1.norm() << std::endl;
    }
  }

  // Check what needs to be updated
  const boost::optional<bool> &update_vbal_cov = params_.update_vbal_cov.value();
  const boost::optional<bool> &update_var = params_.update_var.value();
  const boost::optional<bool> &update_mom = params_.update_mom.value();

  // Load ensemble members sequentially
  if (ensembleConfig1 != boost::none) {
    for (int ie = 0; ie < ens1_ne; ++ie) {
      // Read member
      if (verbosity_) oops::Log::info() <<
      "Info     : -------------------------------------------------------------------" << std::endl;
      if (verbosity_) oops::Log::info() << "Info     : --- Load member " << ie+1 << " / " << ens1_ne
        << std::endl;
      dx1.read(membersConfig1[ie]);

      if (update_vbal_cov != boost::none) {
        if (*update_vbal_cov) {
          // Update vertical covariance
          this->updateVbalCov(dx1.fieldSet(), ie);
        }
      }
      if (update_var != boost::none) {
        if (*update_var) {
          // Update variance
          this->updateVar(dx1.fieldSet(), ie);
        }
      }
      if (update_mom != boost::none) {
        if (*update_mom) {
          // Update moments
          this->updateMom(dx1.fieldSet(), ie, 1);
        }
      }
    }
  }
  if (ensembleConfig2 != boost::none) {
    for (int ie = 0; ie < ens2_ne; ++ie) {
      // Read member
      if (verbosity_) oops::Log::info() <<
      "Info     : -------------------------------------------------------------------" << std::endl;
      if (verbosity_) oops::Log::info() << "Info     : --- Load member " << ie+1 << " / " << ens2_ne
        << std::endl;
      dx2.read(membersConfig2[ie]);
      if (update_mom != boost::none) {
        if (*update_mom) {
          // Update moments
          this->updateMom(dx2.fieldSet(), ie, 2);
        }
      }
    }
  }

  // Run drivers
  this->runDrivers();

  // Partial deallocation
  this->partialDealloc();

  const boost::optional<BUMPOutputNcmpParameters> &outputNcmp = params_.outputNcmp.value();
  const boost::optional<std::vector<BUMPOutputParameters<MODEL>>> &output = params_.output.value();
  if (outputNcmp != boost::none || output != boost::none) {
    // Write parameters
    if (verbosity_) oops::Log::info() <<
    "Info     : -------------------------------------------------------------------" << std::endl;
    if (verbosity_) oops::Log::info() << "Info     : --- Write parameters" << std::endl;
  }
  if (outputNcmp != boost::none) {
    // Open file
    std::ofstream outfile;
    outfile.open(outputNcmp->filepath.value().c_str());

    if (outfile.is_open()) {
      if (verbosity_) oops::Log::test() << "Write number of components in file "
        << outputNcmp->filepath.value() << std::endl;

      // Write parameter
      for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
        for (size_t jvar=0; jvar < activeVarsPerGrid_[jgrid].size(); ++jvar) {
          int ncmp;
          this->getNcmp(jgrid, jvar, ncmp);
          outfile << activeVarsPerGrid_[jgrid][jvar] << ' ' << std::scientific
            << std::setprecision(3) << ncmp << std::endl;
          if (verbosity_) oops::Log::test() << "Number of BUMP output components for "
            << activeVarsPerGrid_[jgrid][jvar] << ": " << ncmp << std::endl;
        }
      }

      // Close file
      outfile.close();
    } else {
      ABORT("BUMP::BUMP: cannot open file");
    }
  }
  if (output != boost::none) {
    for (const auto & outputParam : *output) {
      // Get parameter
      const std::string & param = outputParam.param;

      // Get component
      const int & component = outputParam.component;

      // Select geometry
      if (param == "loc_a_lr"
       || param == "loc_rh_lr"
       || param == "loc_rh1_lr"
       || param == "loc_rh2_lr"
       || param == "loc_rhc_lr"
       || param == "loc_rv_lr"
       || param == "dirac_diag_loc_lr"
       || param == "nicas_norm_lr"
       || param == "dirac_nicas_lr"
       || param == "dirac_nicas_bens_lr") {
        // Get parameter
        dx2.zero(xb.validTime());
        this->getParameter(param, component, 2, dx2.fieldSet());
        dx2.synchronizeFields();

        // Write parameter
        dx2.write(outputParam.incwrite);
        if (verbosity_) oops::Log::test() << "Norm of BUMP output parameter " << param << " - "
          << component << ": " << dx2.norm() << std::endl;
      } else {
        // Get parameter
        dx1.zero(xb.validTime());
        this->getParameter(param, component, 1, dx1.fieldSet());
        dx1.synchronizeFields();

        // Write parameter
        dx1.write(outputParam.incwrite);
        if (verbosity_) oops::Log::test() << "Norm of BUMP output parameter " << param << " - "
          << component << ": "<< dx1.norm() << std::endl;
      }
    }
  }

  // Apply operators
  const boost::optional<std::vector<eckit::LocalConfiguration>>
    &appConfs = params_.appConfs.value();
  if (appConfs != boost::none) {
    if (verbosity_) oops::Log::info() <<
    "Info     : -------------------------------------------------------------------" << std::endl;
    if (verbosity_) oops::Log::info() << "Info     : --- Apply operators" << std::endl;
    if (appConfs->size() > 0) {
      for (const auto & appConf : *appConfs) {
        // Read input file
        eckit::LocalConfiguration inputConf(appConf, "input");
        if (verbosity_) oops::Log::info() << "Info     :        - Input file: " << inputConf
          << std::endl;
        dx1.read(inputConf);

        // Apply BUMP operator
        std::vector<std::string> bumpOperators;
        appConf.get("bump operators", bumpOperators);
        for (const auto & bumpOperator : bumpOperators) {
          if (verbosity_) oops::Log::info() << "Info     :          Apply operator " << bumpOperator
            << std::endl;
          if (bumpOperator == "multiplyVbal") {
            this->multiplyVbal(dx1.fieldSet());
          } else if (bumpOperator == "inverseMultiplyVbal") {
            this->inverseMultiplyVbal(dx1.fieldSet());
          } else if (bumpOperator == "multiplyVbalAd") {
            this->multiplyVbalAd(dx1.fieldSet());
          } else if (bumpOperator == "inverseMultiplyAd") {
            this->inverseMultiplyVbalAd(dx1.fieldSet());
          } else if (bumpOperator == "multiplyStdDev") {
            this->multiplyStdDev(dx1.fieldSet());
          } else if (bumpOperator == "inverseMultiplyStdDev") {
            this->inverseMultiplyStdDev(dx1.fieldSet());
          } else if (bumpOperator == "multiplyNicas") {
            this->multiplyNicas(dx1.fieldSet());
          } else {
              ABORT("BUMP::BUMP: wrong bump operator: " + bumpOperator);
          }
        }

        // ATLAS fieldset to Increment_
        dx1.synchronizeFields();

        // Write file
        eckit::LocalConfiguration outputConf(appConf, "output");
        if (verbosity_) oops::Log::info() << "Info     :          Output file: " << outputConf
          << std::endl;
        dx1.write(outputConf);
      }
    }
  }

  oops::Log::trace() << "BUMP:BUMP constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
BUMP<MODEL>::BUMP(BUMP & other) : keyBUMP_(), activeVarsPerGrid_() {
  for (unsigned int jgrid = 0; jgrid < other.keyBump_.size(); ++jgrid) {
    keyBUMP_.push_back(other.keyBUMP_[jgrid]);
    activeVarsPerGrid_.push_back(other.activeVarsPerGrid_[jgrid]);
  }
  other.keyBUMP_.clear();
  other.activeVarsPerGrid_.clear();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
BUMP<MODEL>::~BUMP() {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    if (keyBUMP_[jgrid] > 0) bump_dealloc_f90(keyBUMP_[jgrid]);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::addMember(const atlas::FieldSet & fset, const int & ie,
                            const int & iens) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_add_member_f90(keyBUMP_[jgrid], fset.get(), ie+1, iens);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::updateVbalCov(const atlas::FieldSet & fset, const int & ie) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_update_vbal_cov_f90(keyBUMP_[jgrid], fset.get(), ie+1);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::updateVar(const atlas::FieldSet & fset, const int & ie) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_update_var_f90(keyBUMP_[jgrid], fset.get(), ie+1);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::updateMom(const atlas::FieldSet & fset, const int & ie,
                            const int & iens) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_update_mom_f90(keyBUMP_[jgrid], fset.get(), ie+1, iens);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::runDrivers() const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_run_drivers_f90(keyBUMP_[jgrid]);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::multiplyVbal(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_vbal_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::inverseMultiplyVbal(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_vbal_inv_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::multiplyVbalAd(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_vbal_ad_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::inverseMultiplyVbalAd(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_vbal_inv_ad_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::multiplyStdDev(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_stddev_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::inverseMultiplyStdDev(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_stddev_inv_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::randomizeNicas(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_randomize_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::multiplyNicas(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_nicas_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::multiplyPsiChiToUV(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_psichi_to_uv_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::multiplyPsiChiToUVAd(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_psichi_to_uv_ad_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::getNcmp(const int & jgrid, const int & jvar, int & ncmp) const {
  bump_get_ncmp_f90(keyBUMP_[jgrid], jvar, ncmp);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::getParameter(const std::string & param, const int & icmp,
  const int & igeom, atlas::FieldSet & fset) const {
  const int npar = param.size();
  const char *cpar = param.c_str();
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_get_parameter_f90(keyBUMP_[jgrid], npar, cpar, icmp, igeom, fset.get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::setNcmp(const int & jgrid, const int & jvar, const int & ncmp) const {
  bump_set_ncmp_f90(keyBUMP_[jgrid], jvar, ncmp);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::setParameter(const std::string & param, const int & icmp,
  const atlas::FieldSet & fset) const {
  const int npar = param.size();
  const char *cpar = param.c_str();
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_set_parameter_f90(keyBUMP_[jgrid], npar, cpar, icmp, fset.get());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP<MODEL>::partialDealloc() const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_partial_dealloc_f90(keyBUMP_[jgrid]);
  }
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_BUMP_BUMP_H_
