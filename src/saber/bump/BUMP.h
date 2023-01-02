/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <math.h>
#include <omp.h>

#include <algorithm>
#include <fstream>
#include <limits>
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

class ValueParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ValueParameters, oops::Parameters)

 public:
  // Variables
  oops::RequiredParameter<std::vector<std::string>> variables{"variables", this};
  // Value
  oops::RequiredParameter<int> value{"value", this};
};

// -----------------------------------------------------------------------------

class GeneralSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeneralSection, oops::Parameters)

 public:
  // Add colors to the log (for display on terminal)
  oops::Parameter<bool> color_log{"color log", false, this};
  // Stream test messages into a dedicated channel
  oops::Parameter<bool> testing{"testing", false, this};
  // Default seed for random numbers
  oops::Parameter<bool> default_seed{"default seed", true, this};
  // Inter-compilers reproducibility
  oops::Parameter<bool> repro_ops{"reproducibility operators", true, this};
  // Reproducibility threshold
  oops::Parameter<double> repro_th{"reproducibility threshold", 1.0e-12, this};
  // Universe radius [in meters]
  oops::Parameter<double> universe_radius{"universe length-scale", 6371229*M_PI, this};
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
  oops::Parameter<std::string> data_directory{"data directory", ".", this};
  // Data prefix
  oops::Parameter<std::string> files_prefix{"files prefix", "", this};
  // Parallel NetCDF I/O
  oops::Parameter<bool> parallel_netcdf{"parallel netcdf", true, this};
  // Number of I/O processors
  oops::Parameter<int> nprocio{"io task number", 20, this};
  // Alias
  oops::Parameter<std::vector<AliasParameter>> alias{"alias", {}, this};
  // Sampling file
  oops::Parameter<std::string> fname_samp{"overriding sampling file", "", this};
  // Vertical covariance files
  oops::Parameter<std::vector<std::string>>
    fname_vbal_cov{"overriding vertical covariance file", {}, this};
  // Vertical balance file
  oops::Parameter<std::string> fname_vbal{"overriding vertical balance file", "", this};
  // Ensemble 1 moments files
  oops::Parameter<std::vector<std::string>> fname_mom{"overriding moments file", {}, this};
  // Ensemble 2 moments files
  oops::Parameter<std::vector<std::string>> fname_mom2{"overriding lowres moments file", {},
    this};
  // NICAS file
  oops::Parameter<std::string> fname_nicas{"overriding nicas file", "", this};
  // Psichitouv transform file
  oops::Parameter<std::string> fname_wind{"overriding psichitouv file", "", this};
};

// -----------------------------------------------------------------------------

class DriversSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DriversSection, oops::Parameters)

 public:
  // Compute covariance, ensemble 1
  oops::Parameter<bool> compute_cov1{"compute covariance", false, this};
  // Compute covariance, ensemble 2
  oops::Parameter<bool> compute_cov2{"compute lowres covariance", false, this};
  // Compute correlation, ensemble 1
  oops::Parameter<bool> compute_cor1{"compute correlation", false, this};
  // Compute correlation, ensemble 2
  oops::Parameter<bool> compute_cor2{"compute lowres correlation", false, this};
  // Compute localization, ensemble 1
  oops::Parameter<bool> compute_loc1{"compute localization", false, this};
  // Compute localization, ensemble 2
  oops::Parameter<bool> compute_loc2{"compute lowres localization", false, this};
  // Compute hybrid weights
  oops::Parameter<bool> compute_hyb{"compute hybrid weights", false, this};
  // Hybrid term source ('randomized static' or 'lowres ensemble')
  oops::Parameter<std::string> hybrid_source{"hybrid source", "", this};
  // Multivariate strategy ('diag_all', 'common', 'common_weighted', 'specific_univariate' or
  // 'specific_multivariate')
  oops::Parameter<std::string> strategy{"multivariate strategy", "", this};
  // Iterative algorithm (ensemble members loaded sequentially)
  oops::Parameter<bool> iterative_algo{"iterative algorithm", false, this};
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
  // Write local vertical covariancee
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
  // Write HDIAG diagnostics
  oops::Parameter<bool> write_hdiag{"write diagnostics", false, this};
  // Write HDIAG components detail
  oops::Parameter<bool> write_hdiag_detail{"write diagnostics detail", false, this};
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
  // Compute wind transform
  oops::Parameter<bool> new_wind{"compute psichitouv", false, this};
  // Read local wind transform
  oops::Parameter<bool> load_wind_local{"read local psichitouv", false, this};
  // Write local wind transform
  oops::Parameter<bool> write_wind_local{"write local psichitouv", false, this};
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
};

// -----------------------------------------------------------------------------

class ModelSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelSection, oops::Parameters)

 public:
  // Level for 2D variables ('first' or 'last')
  oops::Parameter<std::string> lev2d{"level for 2d variables", "first", this};
  // Variables names
  oops::Parameter<std::vector<std::string>> variables{"variables", {}, this};
  // Check that sampling couples and interpolations do not cross mask boundaries
  oops::Parameter<bool> mask_check{"do not cross mask boundaries", false, this};
};

// -----------------------------------------------------------------------------

class EnsembleSizesSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(EnsembleSizesSection, oops::Parameters)

 public:
  // Ensemble 1 size
  oops::Parameter<int> ens1_ne{"total ensemble size", 0, this};
  // Ensemble 1 sub-ensembles number
  oops::Parameter<int> ens1_nsub{"sub-ensembles", 1, this};
  // Ensemble 2 size
  oops::Parameter<int> ens2_ne{"total lowres ensemble size", 0, this};
  // Ensemble 2 sub-ensembles number
  oops::Parameter<int> ens2_nsub{"lowres sub-ensembles", 1, this};
};

// -----------------------------------------------------------------------------

class MaskParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(MaskParameters, oops::Parameters)

 public:
  // Mask restriction type
  oops::RequiredParameter<std::string> mask_type{"type", this};
  // Mask threshold
  oops::Parameter<double> mask_th{"threshold", util::missingValue(double()), this};
  // Mask threshold side ('lower' if mask_th is the lower bound, resp. 'upper')
  oops::Parameter<std::string> mask_lu{"side", "", this};
  // Mask variable
  oops::Parameter<std::string> mask_variable{"variable", "", this};
};

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
  // Activate local diagnostics
  oops::Parameter<bool> local_diag{"local diagnostic", false, this};
  // Local diagnostics calculation radius [in meters]
  oops::Parameter<double> local_rad{"averaging length-scale", 0.0, this};
  // Local diagnostics calculation latitude band half-width [in degrees]
  oops::Parameter<double> local_dlat{"averaging latitude width", 0.0, this};
  // Diagnostic draw type ('random' or 'octahedral')
  oops::Parameter<std::string> draw_type{"grid type", "random", this};
  // Maximum number of random number draws
  oops::Parameter<int> irmax{"max number of draws", 10000, this};
  // Vertical balance C2B to C0A interpolation type ('c0': C0 mesh-based, 'c1': C1 mesh-based
  // or 'si': smooth interpolation)
  oops::Parameter<std::string> interp_type{"interpolation type", "c0", this};
  // Sampling masks
  oops::Parameter<std::vector<MaskParameters>> masks{"masks", {}, this};
  // Threshold on vertically contiguous points for the mask (0 to skip the test)
  oops::Parameter<int> ncontig_th{"contiguous levels threshold", 0, this};
};

// -----------------------------------------------------------------------------

class DiagnosticsSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DiagnosticsSection, oops::Parameters)

 public:
  // Ensemble size
  oops::Parameter<int> ne{"target ensemble size", 0, this};
  // Ensemble size of the hybrid term
  oops::Parameter<int> ne_lr{"target lowres ensemble size", 0, this};
  // Gaussian approximation for asymptotic quantities
  oops::Parameter<bool> gau_approx{"gaussian approximation", false, this};
  // Threshold on generalized kurtosis (3.0 = Gaussian distribution)
  oops::Parameter<double> gen_kurt_th{"generalized kurtosis threshold",
    std::numeric_limits<double>().max(), this};
  // Number of bins for averaged statistics histograms
  oops::Parameter<int> avg_nbins{"histogram bins", 0, this};
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

class VarianceSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VarianceSection, oops::Parameters)

 public:
  // Force specific variance
  oops::Parameter<bool> forced_var{"explicit stddev", false, this};
  // Forced standard-deviation
  oops::Parameter<std::vector<ValueOrProfileParameters>> stddev{"stddev", {}, this};
  // Filter variance
  oops::Parameter<bool> var_filter{"objective filtering", false, this};
  // Number of iterations for the variance filtering (0 for uniform variance)
  oops::Parameter<int> var_niter{"filtering iterations", -1, this};
  // Number of passes for the variance filtering (0 for uniform variance)
  oops::Parameter<int> var_npass{"filtering passes", -1, this};
  // Variance initial filtering support radius [in meters]
  oops::Parameter<std::vector<ValueOrProfileParameters>> var_rhflt{"initial length-scale",
    {}, this};
};

// -----------------------------------------------------------------------------

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

class FitSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FitSection, oops::Parameters)

 public:
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
  oops::Parameter<double> resol{"resolution", 0.0, this};
  // Maximum size of the Sc1 subset
  oops::Parameter<int> nc1max{"max horizontal grid size", 15000, this};
  // NICAS draw type ('random' or 'octahedral')
  oops::Parameter<std::string> nicas_draw_type{"grid type", "random", this};
  // Force specific support radii
  oops::Parameter<bool> forced_radii{"explicit length-scales", false, this};
  // Forced horizontal support radius [in meters]
  oops::Parameter<std::vector<ValueOrProfileParameters>> rh{"horizontal length-scale", {},
    this};
  // Forced vertical support radius
  oops::Parameter<std::vector<ValueOrProfileParameters>> rv{"vertical length-scale", {},
    this};
  // Forced localization weights
  oops::Parameter<std::vector<LocWgtParameters>> loc_wgt{"common localization weights", {},
    this};
  // Minimum level
  oops::Parameter<std::vector<ValueParameters>> min_lev{"minimum level", {}, this};
  // Maximum level
  oops::Parameter<std::vector<ValueParameters>> max_lev{"maximum level", {}, this};
  // NICAS C1B to C0A interpolation type ('c0': C0 mesh-based, 'c1': C1 mesh-based
  // or 'si': smooth interpolation)
  oops::Parameter<std::vector<SpecificTypeParameters>> interp_type{"interpolation type", {},
    this};
  // Positive-definiteness test
  oops::Parameter<bool> pos_def_test{"positive-definiteness test", false, this};
  // Horizontal NICAS interpolation test
  oops::Parameter<bool> interp_test{"horizontal interpolation test", false, this};
};

// -----------------------------------------------------------------------------

class PsichitouvSection : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(PsichitouvSection, oops::Parameters)

 public:
  // Streamfunction variable name
  oops::Parameter<std::string> wind_streamfunction{"stream function", "stream_function", this};
  // Velocity potential variable name
  oops::Parameter<std::string> wind_velocity_potential{"velocity potential", "velocity_potential",
    this};
  // Eastward wind variable name
  oops::Parameter<std::string> wind_eastward{"eastward wind", "eastward_wind", this};
  // Northward wind variable name
  oops::Parameter<std::string> wind_northward{"northward wind", "northward_wind", this};
  // Number of longitudes for the regular grid
  oops::Parameter<int> wind_nlon{"longitudes", 0, this};
  // Number of latitudes for the regular grid
  oops::Parameter<int> wind_nlat{"latitudes", 0, this};
  // Half-width of the Savitzky-Golay to compute derivatives
  oops::Parameter<int> wind_nsg{"savitzky-golay half width", 0, this};
  // Wind inflation to compensate the Savitzky-Golay smoothing
  oops::Parameter<double> wind_inflation{"wind inflation", 1.0, this};
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
  oops::Parameter<GeneralSection> general{"general", GeneralSection(), this};
  // IO parameters
  oops::Parameter<IoSection> io{"io", IoSection(), this};
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
  oops::Parameter<std::vector<LocalProfileParameter>> localProfiles{"local profiles", {}, this};
  // NICAS parameters
  oops::Parameter<NicasSection> nicas{"nicas", NicasSection(), this};
  // Psichitouv parameters
  oops::Parameter<PsichitouvSection> psichitouv{"psichitouv", PsichitouvSection(), this};
  // Dirac parameters
  oops::Parameter<std::vector<DiracPointParameters>> dirac{"dirac", {}, this};

  // External parameters

  // Missing real value
  oops::Parameter<double> msvalr{"msvalr", util::missingValue(double()), this};
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
