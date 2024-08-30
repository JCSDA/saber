/*
 * (C) Crown Copyright 2022-2024 Met Office
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/spectralb/SpectralCovariance.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"

#include "eckit/exception/Exceptions.h"

#include "oops/mpi/mpi.h"
#include "oops/util/AtlasArrayUtil.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

#include "saber/oops/Utilities.h"
#include "saber/spectralb/CovarianceStatisticsUtils.h"
#include "saber/util/BinnedFieldSetHelpers.h"
#include "saber/util/Calibration.h"

// -----------------------------------------------------------------------------

namespace saber {
namespace spectralb {

namespace {

eckit::LocalConfiguration createNetCDFHeaderInput(
    const SpectralCovarianceParameters & params,
    const oops::Variables & activeVars) {
  const std::string statsType = "spectral vertical covariance";
  const std::string binType = "total wavenumber";
  bool doingCalibration(params.calibrationParams.value() != boost::none);
  const eckit::LocalConfiguration conf = params.toConfiguration();
  eckit::LocalConfiguration netCDFConf;

  util::createCalibrationNetCDFHeaderInput(conf, statsType, binType,
                                           activeVars, doingCalibration,
                                           netCDFConf);
  return netCDFConf;
}

}  // namespace

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<SpectralCovariance> makerSpectralCovariance_("spectral covariance");

// -----------------------------------------------------------------------------
SpectralCovariance::SpectralCovariance(const oops::GeometryData & geometryData,
                                       const oops::Variables & centralVars,
                                       const eckit::Configuration & covarConf,
                                       const Parameters_ & params,
                                       const oops::FieldSet3D & xb,
                                       const oops::FieldSet3D & fg)
  : SaberCentralBlockBase(params, xb.validTime()), params_(params),
    activeVars_(getActiveVars(params, centralVars)),
    netCDFConf_(createNetCDFHeaderInput(params, activeVars_)),
    spectralVerticalCovariances_(),
    geometryData_(geometryData),
    specFunctionSpace_(geometryData_.functionSpace()),
    datetime_(xb.validTime()),
    sample_size_(0)
{
  oops::Log::trace() << classname() << "::SpectralCovariance starting " << std::endl;

  oops::Log::trace() << classname() << "::SpectralCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralCovariance::randomize(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;

  oops::Log::error() << "randomization with spectral covariance saber block"
                     << " is not supported. Instead please use 'ID' central block"
                     << " and 'square root of spectral covariance' outer block."
                     << std::endl;
  throw(eckit::FunctionalityNotSupported(
        "use ID and square root of spectral covariance instead.", Here()));
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralCovariance::multiply(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  if (params_.skipVerticalConv.value()) {
    specutils::spectralHorizontalFilter(activeVars_,
                                        specFunctionSpace_,
                                        spectralVerticalCovariances_,
                                        fieldSet.fieldSet());
  } else {
    specutils::spectralVerticalConvolution(activeVars_,
                                           specFunctionSpace_,
                                           spectralVerticalCovariances_,
                                           fieldSet.fieldSet());
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralCovariance::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  // Note that the read can occur in either calibration mode
  // using calibrationReadParams or in standard covariance mode using
  // readParams
  spectralbReadParameters sparams;
  const auto & calibparams = params_.calibrationParams.value();
  if (calibparams != boost::none) {
    const auto & calibrationReadParams = calibparams->calibrationReadParams.value();
    if (calibrationReadParams != boost::none) {
       sparams = calibrationReadParams.value();
    }
  } else {
    sparams = *params_.readParams.value();
  }

  const auto & umatrixNetCDFParams = sparams.umatrixNetCDFNames.value();

  const int nSpectralBins = specFunctionSpace_.truncation() + 1;  // 2N
  spectralVerticalCovariances_.clear();

  for (std::size_t i = 0; i < activeVars_.size(); ++i) {
    //  allocate vert cov field based on activeVars and spectralfunctionspace_
    auto spectralVertCov =
      atlas::Field(activeVars_[i].name(),
                   atlas::array::make_datatype<double>(),
                   atlas::array::make_shape(nSpectralBins,
                                            activeVars_[i].getLevels(),
                                            activeVars_[i].getLevels()));
    if (umatrixNetCDFParams != boost::none) {
      const oops::Variables netCDFVars(umatrixNetCDFParams.value());
      specutils::createSpectralCovarianceFromUMatrixFile(activeVars_[i].name(),
                                                         netCDFVars[i].name(),
                                                         sparams,
                                                         spectralVertCov);
    } else {
      specutils::readSpectralCovarianceFromFile(activeVars_[i].name(),
                                                activeVars_[i].name(),
                                                sparams,
                                                spectralVertCov);
    }

    spectralVerticalCovariances_.add(spectralVertCov);
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}


void SpectralCovariance::directCalibration(const oops::FieldSets &
                                           MOSpectralCovariancesEns) {
  oops::Log::trace() << classname() << "::directCalibration starting" << std::endl;

  const auto & calibparams = params_.calibrationParams.value();
  ASSERT(calibparams != boost::none);
  const auto & calibrationReadParams = calibparams->calibrationReadParams.value();
  if (calibrationReadParams != boost::none) {
    SpectralCovariance::read();
  } else {
    oops::Log::info() << "not reading cov file but allocating and zeroing" << std::endl;
    spectralVerticalCovariances_.clear();
    const int nSpectralBins = specFunctionSpace_.truncation() + 1;  // 2N
    for (std::size_t i = 0; i < activeVars_.size(); ++i) {
      auto spectralVertCov = atlas::Field(activeVars_[i].name(),
                                          atlas::array::make_datatype<double>(),
        atlas::array::make_shape(nSpectralBins,
                                 activeVars_[i].getLevels(),
                                 activeVars_[i].getLevels()));

      auto spectralVertCovView = atlas::array::make_view<double, 3>(spectralVertCov);
      spectralVertCovView.assign(0.0);
      spectralVerticalCovariances_.add(spectralVertCov);
    }
  }

  if (MOSpectralCovariancesEns.ens_size() > 0) {
    sample_size_ = specutils::updateSpectralVerticalCovariances(
      MOSpectralCovariancesEns, sample_size_, spectralVerticalCovariances_);
  }

  oops::Log::trace() << classname() << "::directCalibration done" << std::endl;
}

void SpectralCovariance::write() const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  util::calibrationWriteParameters writeParams =
    (params_.calibrationParams.value())->writeParams.value();
  eckit::LocalConfiguration writeConfig;
  writeParams.serialize(writeConfig);

  const std::string filepath{writeParams.filePath};
  const std::string mpi_pattern = writeParams.mpiPattern;
  const std::string mpi_size = std::to_string(eckit::mpi::comm().size());
  ::util::seekAndReplace(writeConfig, mpi_pattern, mpi_size);
  const std::string ncfilepath = "./" + writeConfig.getString("file path");

  using atlas::array::make_view;
  using atlas::idx_t;

  // The spectralVerticalCovariances that we write should be a gathered version of
  // the one in memory  ... it should not affect the internal version.
  // gather and sum on pe 0
  const std::size_t root(0);
  atlas::FieldSet spectralVertCovToWrite = util::gatherSumFieldSet(geometryData_.comm(),
                                                                   root,
                                                                   spectralVerticalCovariances_);
  const std::vector<std::string> dim_names{"binning_index",
                                           "levels_index_1",
                                           "levels_index_2"};
  const std::vector<atlas::idx_t> dim_sizes{spectralVertCovToWrite[0].shape()[0],
                                            spectralVertCovToWrite[0].shape()[1],
                                            spectralVertCovToWrite[0].shape()[2]};
  oops::Variables fsetVars;
  std::vector<std::vector<std::string>> dim_names_for_every_var;
  for (const std::string & name : netCDFConf_.keys()) {
    const auto varname1 =
      util::getAttributeValue<std::string>(netCDFConf_, name, "variable_name_1");
    const auto varname2 =
      util::getAttributeValue<std::string>(netCDFConf_, name, "variable_name_2");
    ASSERT(varname1 == varname2);
    fsetVars.push_back(oops::Variable{name, oops::VariableMetaData(),
                                      activeVars_[varname1].getLevels()});
    spectralVertCovToWrite[varname1].rename(name);
    dim_names_for_every_var.push_back(dim_names);
  }

  eckit::LocalConfiguration netCDFConf = netCDFConf_;
  util::setAttribute<std::string>(
    netCDFConf, "global metadata", "covariance_name", "string", writeParams.covName);
  util::setAttribute<std::string>(
    netCDFConf, "global metadata", "date_time", "string", datetime_.toString());
  util::setAttribute<std::int32_t>(
    netCDFConf, "global metadata", "no_of_samples", "int32",
    static_cast<std::int32_t>(sample_size_));

  std::vector<int> netcdf_general_ids;
  std::vector<int> netcdf_dim_ids;
  std::vector<int> netcdf_var_ids;
  std::vector<std::vector<int>> netcdf_dim_var_ids;

  if (eckit::mpi::comm().rank() == root) {
    ::util::atlasArrayWriteHeader(ncfilepath,
                                  dim_names,
                                  dim_sizes,
                                  fsetVars,
                                  dim_names_for_every_var,
                                  netCDFConf,
                                  netcdf_general_ids,
                                  netcdf_dim_ids,
                                  netcdf_var_ids,
                                  netcdf_dim_var_ids);

    std::size_t t(0);
    for (const atlas::Field & fld : spectralVertCovToWrite) {
      auto fview = atlas::array::make_view<const double, 3>(fld);
      ::util::atlasArrayWriteData(netcdf_general_ids,
                                  netcdf_var_ids[t],
                                  fview);
      ++t;
    }
  }

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralCovariance::print(std::ostream & os) const {
  os << classname();
}

}  // namespace spectralb
}  // namespace saber
