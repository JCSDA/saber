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
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/bump/type_bump.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

class ValueOrProfileParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ValueOrProfileParameters, oops::Parameters)

 public:
  // Variables
  oops::RequiredParameter<std::vector<std::string>> variables{"variables", this};
  // Value
  oops::OptionalParameter<double> value{"value", this};
  // Profile
  oops::OptionalParameter<std::vector<double>> profile{"profile", this};
};

// -----------------------------------------------------------------------------

class GeneralSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeneralSection, oops::Parameters)

 public:
  // Add colors to the log (for display on terminal)
  oops::OptionalParameter<bool> color_log{"color log", this};
  // Stream test messages into a dedicated channel
  oops::OptionalParameter<bool> testing{"testing", this};
  // Default seed for random numbers
  oops::OptionalParameter<bool> default_seed{"default seed", this};
  // Inter-compilers reproducibility
  oops::OptionalParameter<bool> repro_ops{"reproducibility operators", this};
  // Reproducibility threshold
  oops::OptionalParameter<double> repro_th{"reproducibility threshold", this};
  // Universe radius [in meters]
  oops::OptionalParameter<double> universe_radius{"universe length-scale", this};
};

// -----------------------------------------------------------------------------

class AliasParameter : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(AliasParameter, oops::Parameters)

 public:
  // In code
  oops::RequiredParameter<std::string> in_code{"in code", this};
  // In file
  oops::RequiredParameter<std::string> in_file{"in file", this};
};

class IoSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(IoSection, oops::Parameters)

 public:
  // Data directory
  oops::OptionalParameter<std::string> data_directory{"data directory", this};
  // Data prefix
  oops::OptionalParameter<std::string> files_prefix{"files prefix", this};
  // Parallel NetCDF I/O
  oops::OptionalParameter<bool> parallel_netcdf{"parallel netcdf", this};
  // Number of I/O processors
  oops::OptionalParameter<int> nprocio{"io task number", this};
  // Alias
  oops::OptionalParameter<std::vector<AliasParameter>> alias{"alias", this};
  // Sampling file
  oops::OptionalParameter<std::string> fname_samp{"overriding sampling file", this};
  // Vertical covariance files
  oops::OptionalParameter<std::vector<std::string>>
    fname_vbal_cov{"overriding vertical covariance file", this};
  // Vertical balance file
  oops::OptionalParameter<std::string> fname_vbal{"overriding vertical balance file", this};
  // Moments files
  oops::OptionalParameter<std::vector<std::string>> fname_mom{"overriding moments file", this};
  // NICAS file
  oops::OptionalParameter<std::string> fname_nicas{"overriding nicas file", this};
  // Wind transform file
  oops::OptionalParameter<std::string> fname_wind{"overriding psichitouv file", this};
};

// -----------------------------------------------------------------------------

class DriversSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DriversSection, oops::Parameters)

 public:
  // Compute covariance, ensemble 1
  oops::OptionalParameter<bool> compute_cov1{"compute covariance", this};
  // Compute covariance, ensemble 2
  oops::OptionalParameter<bool> compute_cov2{"compute lowres covariance", this};
  // Compute correlation, ensemble 1
  oops::OptionalParameter<bool> compute_cor1{"compute correlation", this};
  // Compute correlation, ensemble 2
  oops::OptionalParameter<bool> compute_cor2{"compute lowres correlation", this};
  // Compute localization, ensemble 1
  oops::OptionalParameter<bool> compute_loc1{"compute localization", this};
  // Compute localization, ensemble 2
  oops::OptionalParameter<bool> compute_loc2{"compute lowres localization", this};
  // Compute hybrid weights
  oops::OptionalParameter<bool> compute_hyb{"compute hybrid weights", this};
  // Hybrid term source ('randomized static' or 'lowres ensemble')
  oops::OptionalParameter<std::string> hybrid_source{"hybrid source", this};
  // Multivariate strategy ('diag_all', 'common', 'common_weighted', 'specific_univariate' or
  // 'specific_multivariate')
  oops::OptionalParameter<std::string> strategy{"multivariate strategy", this};
  // Iterative algorithm (ensemble members loaded sequentially)
  oops::OptionalParameter<bool> iterative_algo{"iterative algorithm", this};
  // New normality test
  oops::OptionalParameter<bool> new_normality{"compute normality", this};
  // Load local sampling
  oops::OptionalParameter<bool> load_samp_local{"load local sampling", this};
  // Load global sampling
  oops::OptionalParameter<bool> load_samp_global{"load global sampling", this};
  // Write local sampling
  oops::OptionalParameter<bool> write_samp_local{"write local sampling", this};
  // Write global sampling
  oops::OptionalParameter<bool> write_samp_global{"write global sampling", this};
  // Write sampling grids
  oops::OptionalParameter<bool> write_samp_grids{"write sampling grids", this};
  // New vertical covariance
  oops::OptionalParameter<bool> new_vbal_cov{"compute vertical covariance", this};
  // Load local vertical covariance
  oops::OptionalParameter<bool> load_vbal_cov{"load vertical covariance", this};
  // Write local vertical covariancee
  oops::OptionalParameter<bool> write_vbal_cov{"write vertical covariance", this};
  // Compute vertical balance operator
  oops::OptionalParameter<bool> new_vbal{"compute vertical balance", this};
  // Load local vertical balance operator
  oops::OptionalParameter<bool> load_vbal{"load vertical balance", this};
  // Write vertical balance operator
  oops::OptionalParameter<bool> write_vbal{"write vertical balance", this};
  // Compute variance
  oops::OptionalParameter<bool> new_var{"compute variance", this};
  // Compute moments
  oops::OptionalParameter<bool> new_mom{"compute moments", this};
  // Load sampling moments
  oops::OptionalParameter<bool> load_mom{"load moments", this};
  // Write sampling moments
  oops::OptionalParameter<bool> write_mom{"write moments", this};
  // Compute HDIAG
  oops::OptionalParameter<bool> new_hdiag{"compute diagnostics", this};
  // Write HDIAG diagnostics
  oops::OptionalParameter<bool> write_hdiag{"write diagnostics", this};
  // Write HDIAG components detail
  oops::OptionalParameter<bool> write_hdiag_detail{"write diagnostics detail", this};
  // Compute NICAS
  oops::OptionalParameter<bool> new_nicas{"compute nicas", this};
  // Load local NICAS parameters
  oops::OptionalParameter<bool> load_nicas_local{"load local nicas", this};
  // Load global NICAS parameters
  oops::OptionalParameter<bool> load_nicas_global{"load global nicas", this};
  // Write local NICAS parameters
  oops::OptionalParameter<bool> write_nicas_local{"write local nicas", this};
  // Write global NICAS parameters
  oops::OptionalParameter<bool> write_nicas_global{"write global nicas", this};
  // Write NICAS grids
  oops::OptionalParameter<bool> write_nicas_grids{"write nicas grids", this};
  // Compute wind transform
  oops::OptionalParameter<bool> new_wind{"compute psichitouv", this};
  // Load local wind transform
  oops::OptionalParameter<bool> load_wind_local{"load local psichitouv", this};
  // Write local wind transform
  oops::OptionalParameter<bool> write_wind_local{"write local psichitouv", this};
  // Test vertical balance inverse
  oops::OptionalParameter<bool> check_vbal{"vertical balance inverse test", this};
  // Test adjoints
  oops::OptionalParameter<bool> check_adjoints{"adjoints test", this};
  // Test NICAS normalization (number of tests)
  oops::OptionalParameter<int> check_normalization{"normalization test", this};
  // Test NICAS application on diracs
  oops::OptionalParameter<bool> check_dirac{"internal dirac test", this};
  // Test NICAS randomization
  oops::OptionalParameter<bool> check_randomization{"randomization test", this};
  // Test HDIAG-NICAS consistency
  oops::OptionalParameter<bool> check_consistency{"internal consistency test", this};
  // Test HDIAG optimality
  oops::OptionalParameter<bool> check_optimality{"localization optimality test", this};
};

// -----------------------------------------------------------------------------

class ModelSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelSection, oops::Parameters)

 public:
  // Level for 2D variables ('first' or 'last')
  oops::OptionalParameter<std::string> lev2d{"level for 2d variables", this};
  // Variables names
  oops::OptionalParameter<std::vector<std::string>> variables{"variables", this};
  // Check that sampling couples and interpolations do not cross mask boundaries
  oops::OptionalParameter<bool> mask_check{"do not cross mask boundaries", this};
};

// -----------------------------------------------------------------------------

class EnsembleSizesSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(EnsembleSizesSection, oops::Parameters)

 public:
  // Ensemble 1 size
  oops::OptionalParameter<int> ens1_ne{"total ensemble size", this};
  // Ensemble 1 sub-ensembles number
  oops::OptionalParameter<int> ens1_nsub{"sub-ensembles", this};
  // Ensemble 2 size
  oops::OptionalParameter<int> ens2_ne{"total lowres ensemble size", this};
  // Ensemble 2 sub-ensembles number
  oops::OptionalParameter<int> ens2_nsub{"lowres sub-ensembles", this};
};

// -----------------------------------------------------------------------------

class MaskParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(MaskParameters, oops::Parameters)

 public:
  // Mask restriction type
  oops::RequiredParameter<std::string> mask_type{"type", this};
  // Mask threshold
  oops::OptionalParameter<double> mask_th{"threshold", this};
  // Mask threshold side ('lower' if mask_th is the lower bound, resp. 'upper')
  oops::OptionalParameter<std::string> mask_lu{"side", this};
  // Mask variable
  oops::OptionalParameter<std::string> mask_variable{"variable", this};
};

class SamplingSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SamplingSection, oops::Parameters)

 public:
  // Computation grid size
  oops::OptionalParameter<int> nc1{"computation grid size", this};
  // Diagnostic grid size
  oops::OptionalParameter<int> nc2{"diagnostic grid size", this};
  // Number of distance classes
  oops::OptionalParameter<int> nc3{"distance classes", this};
  // Number of angular sectors
  oops::OptionalParameter<int> nc4{"angular sectors", this};
  // Class size (for sam_type='hor'), should be larger than the typical grid cell size [in meters]
  oops::OptionalParameter<double> dc{"distance class width", this};
  // Reduced number of levels for diagnostics
  oops::OptionalParameter<int> nl0r{"reduced levels", this};
  // Activate local diagnostics
  oops::OptionalParameter<bool> local_diag{"local diagnostic", this};
  // Local diagnostics calculation radius [in meters]
  oops::OptionalParameter<double> local_rad{"averaging length-scale", this};
  // Local diagnostics calculation latitude band half-width [in degrees]
  oops::OptionalParameter<double> local_dlat{"averaging latitude width", this};
  // Diagnostic draw type ('random' or 'octahedral')
  oops::OptionalParameter<std::string> draw_type{"grid type", this};
  // Maximum number of random number draws
  oops::OptionalParameter<int> irmax{"max number of draws", this};
  // Vertical balance C2B to C0A interpolation type ('c0': C0 mesh-based, 'c1': C1 mesh-based
  // or 'si': smooth interpolation)
  oops::OptionalParameter<std::string> interp_type{"interpolation type", this};
  // Sampling masks
  oops::OptionalParameter<std::vector<MaskParameters>> masks{"masks", this};
  // Threshold on vertically contiguous points for the mask (0 to skip the test)
  oops::OptionalParameter<int> ncontig_th{"contiguous levels threshold", this};
};

// -----------------------------------------------------------------------------

class LocalizationSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LocalizationSection, oops::Parameters)

 public:
  // Ensemble size
  oops::OptionalParameter<int> ne{"target ensemble size", this};
  // Ensemble size of the hybrid term
  oops::OptionalParameter<int> ne_lr{"target lowres ensemble size", this};
  // Gaussian approximation for asymptotic quantities
  oops::OptionalParameter<bool> gau_approx{"gaussian approximation", this};
  // Threshold on generalized kurtosis (3.0 = Gaussian distribution)
  oops::OptionalParameter<double> gen_kurt_th{"generalized kurtosis threshold", this};
  // Number of bins for averaged statistics histograms
  oops::OptionalParameter<int> avg_nbins{"histogram bins", this};
};

// -----------------------------------------------------------------------------

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

class VerticalBalanceSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VerticalBalanceSection, oops::Parameters)

 public:
  // Vertical balance parameters
  oops::OptionalParameter<std::vector<VerticalBalanceBlockParameters>> vbal{"vbal", this};
  // Pseudo-inverse for auto-covariance
  oops::OptionalParameter<bool> vbal_pseudo_inv{"pseudo inverse", this};
  // Dominant mode for pseudo-inverse
  oops::OptionalParameter<int> vbal_pseudo_inv_mmax{"dominant mode", this};
  // Variance threshold to compute the dominant mode for pseudo-inverse
  oops::OptionalParameter<double> vbal_pseudo_inv_var_th{"variance threshold", this};
  // Identity vertical balance for tests
  oops::OptionalParameter<bool> vbal_id{"identity blocks", this};
};

// -----------------------------------------------------------------------------

class VarianceSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VarianceSection, oops::Parameters)

 public:
  // Force specific variance
  oops::OptionalParameter<bool> forced_var{"explicit stddev", this};
  // Forced standard-deviation
  oops::OptionalParameter<std::vector<ValueOrProfileParameters>> stddev{"stddev", this};
  // Filter variance
  oops::OptionalParameter<bool> var_filter{"objective filtering", this};
  // Number of iterations for the variance filtering (0 for uniform variance)
  oops::OptionalParameter<int> var_niter{"filtering iterations", this};
  // Number of passes for the variance filtering (0 for uniform variance)
  oops::OptionalParameter<int> var_npass{"filtering passes", this};
  // Variance initial filtering support radius [in meters]
  oops::OptionalParameter<std::vector<ValueOrProfileParameters>> var_rhflt{"initial length-scale",
    this};
};

// -----------------------------------------------------------------------------

class OptimalityTestSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OptimalityTestSection, oops::Parameters)

 public:
  // Number of length-scale factors for optimization
  oops::OptionalParameter<int> optimality_nfac{"half number of factors", this};
  // Increments of length-scale factors for optimization
  oops::OptionalParameter<double> optimality_delta{"factors increment", this};
  // Number of test vectors for optimization
  oops::OptionalParameter<int> optimality_ntest{"test vectors", this};
};

// -----------------------------------------------------------------------------

class FitSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FitSection, oops::Parameters)

 public:
  // Horizontal filtering suport radius [in meters]
  oops::OptionalParameter<double> diag_rhflt{"horizontal filtering length-scale", this};
  // Vertical filtering support radius
  oops::OptionalParameter<double> diag_rvflt{"vertical filtering length-scale", this};
  // Number of levels between interpolation levels
  oops::OptionalParameter<int> fit_dl0{"vertical stride", this};
  // Number of components in the fit function
  oops::OptionalParameter<int> fit_ncmp{"number of components", this};
};

// -----------------------------------------------------------------------------

class LocalProfileParameter : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LocalProfileParameter, oops::Parameters)

 public:
  // Longitudes of the local diagnostics profiles to write [in degrees]
  oops::RequiredParameter<double> lon_ldwv{"longitude", this};
  // Latitudes of the local diagnostics profiles to write [in degrees]
  oops::RequiredParameter<double> lat_ldwv{"latitude", this};
  // Name of the local diagnostics profiles to write
  oops::RequiredParameter<std::string> name_ldwv{"name", this};
};

// -----------------------------------------------------------------------------

class SpecificTypeParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SpecificTypeParameters, oops::Parameters)

 public:
  // Resolution
  oops::RequiredParameter<std::vector<std::string>> variables{"variables", this};
  // Type
  oops::RequiredParameter<std::string> type{"type", this};
};

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

class NicasSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(NicasSection, oops::Parameters)

 public:
  // Resolution
  oops::OptionalParameter<double> resol{"resol", this};
  // Maximum size of the Sc1 subset
  oops::OptionalParameter<int> nc1max{"nc1max", this};
  // NICAS draw type ('random' or 'octahedral')
  oops::OptionalParameter<std::string> nicas_draw_type{"nicas_draw_type", this};
  // Force specific support radii
  oops::OptionalParameter<bool> forced_radii{"forced_radii", this};
  // Forced horizontal support radius [in meters]
  oops::OptionalParameter<std::vector<ValueOrProfileParameters>> rh{"rh", this};
  // Forced vertical support radius
  oops::OptionalParameter<std::vector<ValueOrProfileParameters>> rv{"rv", this};
  // Forced localization weights
  oops::OptionalParameter<std::vector<LocWgtParameters>> loc_wgt{"loc_wgt", this};
  // Minimum level
  oops::OptionalParameter<std::vector<ValueOrProfileParameters>> min_lev{"min_lev", this};
  // Maximum level
  oops::OptionalParameter<std::vector<ValueOrProfileParameters>> max_lev{"max_lev", this};
  // NICAS C1B to C0A interpolation type ('c0': C0 mesh-based, 'c1': C1 mesh-based
  // or 'si': smooth interpolation)
  oops::OptionalParameter<std::vector<SpecificTypeParameters>> interp_type{"interp_type", this};
  // Positive-definiteness test
  oops::OptionalParameter<bool> pos_def_test{"pos_def_test", this};
  // Horizontal NICAS interpolation test
  oops::OptionalParameter<bool> interp_test{"interp_test", this};
};

// -----------------------------------------------------------------------------

class WindSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(WindSection, oops::Parameters)

 public:
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

class BUMPParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BUMPParameters, oops::Parameters)

 public:
  // Internal parameters

  // General parameters
  oops::OptionalParameter<GeneralSection> general{"general", this};
  // Drivers parameters
  oops::OptionalParameter<DriversSection> drivers{"drivers", this};
  // IO parameters
  oops::OptionalParameter<IoSection> io{"io", this};
  // Model parameters
  oops::OptionalParameter<ModelSection> model{"model", this};
  // Ensemble sizes parameters
  oops::OptionalParameter<EnsembleSizesSection> ensembleSizes{"ensemble sizes", this};
  // Sampling parameters
  oops::OptionalParameter<SamplingSection> sampling{"sampling", this};
  // Localization parameters
  oops::OptionalParameter<LocalizationSection> localization{"localization", this};
  // Vertical balance parameters
  oops::OptionalParameter<VerticalBalanceSection> verticalBalance{"vertical balance", this};
  // Variance parameters
  oops::OptionalParameter<VarianceSection> variance{"variance", this};
  // Optimality test parameters
  oops::OptionalParameter<OptimalityTestSection> optimalityTest{"optimality test", this};
  // Fit parameters
  oops::OptionalParameter<FitSection> fit{"fit", this};
  // Local profiles parameters
  oops::OptionalParameter<std::vector<LocalProfileParameter>> localProfiles{"local profiles",
    this};
  // NICAS parameters
  oops::OptionalParameter<NicasSection> nicas{"nicas", this};
  // Wind parameters
  oops::OptionalParameter<WindSection> wind{"wind", this};
  // Dirac parameters
  oops::OptionalParameter<std::vector<DiracPointParameters>> dirac{"dirac", this};

  // External parameters

  // Ensemble 1 parameters
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble1{"ensemble", this};
  // Ensemble 2 parameters
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble2{"lowres ensemble", this};
  // Grids
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> grids{"grids", this};
  // Operators application
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> appConfs{"operators application",
    this};
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
       const size_t & ens1_ne_in = 0,
       const atlas::FunctionSpace & functionSpace2 = NULL,
       const atlas::FieldSet & extraFields2 = NULL,
       const std::vector<atlas::FieldSet> & fsetVec2 = {},
       const size_t & ens2_ne_in = 0);

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
  void finalize() const;

 private:
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
