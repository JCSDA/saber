/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

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

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

#include "saber/bump/type_bump.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

class BUMPParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BUMPParameters, oops::Parameters)

 public:
  // External parameters

  // Ensemble 1 parameters
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble1{"ensemble", this};
  // Ensemble 2 parameters
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble2{"lowres ensemble", this};
  // Missing value (real)
  oops::OptionalParameter<double> msvalr{"msvalr", this};
  // Grids
  oops::OptionalParameter<eckit::LocalConfiguration> grids{"grids", this};
  // Operators application
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> appConfs{"operators application",
    this};

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
};

// -----------------------------------------------------------------------------

class BUMP {
 public:
  // Constructor
  BUMP(const eckit::mpi::Comm &,
       const atlas::FunctionSpace &,
       const atlas::FieldSet &,
       const std::vector<size_t> &,
       const oops::Variables &,
       const BUMPParameters &,
       const std::vector<atlas::FieldSet> &,
       const atlas::FunctionSpace & functionSpace2 = NULL,
       const atlas::FieldSet & extraFields2 = NULL,
       const std::vector<atlas::FieldSet> & fsetVec2 = {},
       const size_t & ens1_ne_in = 0,
       const size_t & ens2_ne_in = 0);

  // Copy-constructor
  explicit BUMP(BUMP &);

  // Destructor
  ~BUMP();

  // Accessors
  const std::vector<eckit::LocalConfiguration> memberConfig1() const {return membersConfig1_;}
  const std::vector<eckit::LocalConfiguration> memberConfig2() const {return membersConfig2_;}

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
  void getParameter(const std::string &, const int &, const int &, atlas::FieldSet &) const;
  void setNcmp(const int &, const int &) const;
  void setParameter(const std::string &, const int &, const atlas::FieldSet &) const;
  void partialDealloc() const;

 private:
  bool verbosity_;
  BUMPParameters params_;
  const oops::Variables activeVars_;
  std::vector<int> keyBUMP_;
  std::vector<eckit::LocalConfiguration> membersConfig1_;
  std::vector<eckit::LocalConfiguration> membersConfig2_;
  std::vector<oops::Variables> activeVarsPerGrid_;
};

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
