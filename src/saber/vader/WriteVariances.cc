/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/WriteVariances.h"

#include <netcdf.h>
#include <sys/stat.h>

#include <sstream>
#include <utility>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/util/Earth.h"
#include "atlas/util/GaussianLatitudes.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/FieldSet3D.h"
#include "oops/util/AtlasArrayUtil.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

#include "saber/util/BinnedFieldSetHelpers.h"
#include "saber/util/Calibration.h"
#include "saber/vader/VarianceAccumulationUtils.h"

namespace saber {
namespace vader {

namespace  {

std::size_t getSizeOwned(const atlas::FunctionSpace & fs) {
  std::size_t sizeOwned = 0;
  bool fsViable(false);
  if (fs.type() == "StructuredColumns") {
    sizeOwned = atlas::functionspace::StructuredColumns(fs).sizeOwned();
    fsViable = true;
  } else if (fs.type() == "NodeColumns") {
    std::string s = atlas::functionspace::NodeColumns(fs).mesh().grid().name();
    if (s.substr(0, 6) == "CS-LFR") {
      sizeOwned = atlas::functionspace::CubedSphereNodeColumns(fs).sizeOwned();
      fsViable = true;
    }
  }
  ASSERT(fsViable);
  return sizeOwned;
}

std::string getBinType(const WriteVariancesParameters & params) {
  binningParameters bparams = params.binning;
  const std::string binType = bparams.type;
  return binType;
}


eckit::LocalConfiguration createNetCDFHeaderInput(const std::string & binType,
                                                  const WriteVariancesParameters & params,
                                                  const oops::Variables & vars) {
  const std::string statsType = params.statisticsType;
  const std::vector<std::string> fnames(params.fieldNames);
  const bool doingCalibration(params.calibrationParams.value() != boost::none);
  oops::Variables fvars;
  for (const std::string & fname : fnames) {
    fvars.push_back(vars[fname]);
  }

  const eckit::LocalConfiguration conf = params.toConfiguration();
  eckit::LocalConfiguration netCDFConf;

  // global header is deferred to "direct calibration" to get the number of ensemble members.
  util::setAttribute<std::string>(netCDFConf, "longitude", "long_name",
                                  "string", "longitude");
  util::setAttribute<std::string>(netCDFConf, "longitude", "units",
                                  "string", "degrees");
  util::setAttribute<std::string>(netCDFConf, "latitude", "long_name",
                                  "string", "latitude");
  util::setAttribute<std::string>(netCDFConf, "latitude", "units",
                                  "string", "degrees");
  std::string binWeights = binType + " weights";
  util::setAttribute<std::string>(netCDFConf, binWeights, "long_name",
                                  "string", binWeights);
  util::setAttribute<std::string>(netCDFConf, binWeights, "binning type",
                                  "string", binType);
  std::string binIndices = binType + " indices";
  util::setAttribute<std::string>(netCDFConf, binIndices, "long_name",
                                  "string", binIndices);
  util::setAttribute<std::string>(netCDFConf, binIndices, "binning type",
                                  "string", binType);

  util::createCalibrationNetCDFHeaderInput(conf, statsType, binType, fvars,
                                           doingCalibration, netCDFConf);

  return netCDFConf;
}

atlas::FieldSet createHorizontalGridPointWeights(const std::string & binType,
                                                 const std::size_t sizeOwned,
                                                 const atlas::FunctionSpace & fs) {
  atlas::Field gridAreaWgt =
    atlas::Field(binType + " weights",  atlas::array::make_datatype<double>(),
                           atlas::array::make_shape(sizeOwned, 1));
  const std::string binTypeIdx = binType + " local indices";
  atlas::Field gridAreaIdx =
    atlas::Field(binTypeIdx, atlas::array::make_datatype<std::int32_t>(),
                 atlas::array::make_shape(sizeOwned, 1));

  const std::string binTypeGlIdx = binType + " global bins";
  atlas::Field binGlobalIdx =
    atlas::Field(binTypeGlIdx, atlas::array::make_datatype<std::int32_t>(),
                 atlas::array::make_shape(sizeOwned));

  auto gridAreaWgtView = atlas::array::make_view<double, 2>(gridAreaWgt);
  auto gridAreaIdxView = atlas::array::make_view<std::int32_t, 2>(gridAreaIdx);
  auto globalIdxView = atlas::array::make_view<atlas::gidx_t, 1>(fs.global_index());
  auto binGlobalIdxView = atlas::array::make_view<std::int32_t, 1>(binGlobalIdx);

  for (std::size_t jn = 0; jn < sizeOwned; jn++) {
    gridAreaWgtView(jn, 0) = 1.0;
    gridAreaIdxView(jn, 0) = static_cast<std::int32_t>(jn);
    binGlobalIdxView(jn) = static_cast<std::int32_t>(globalIdxView(jn)-1);
  }
  atlas::FieldSet binningData;
  binningData.add(gridAreaWgt);
  binningData.add(gridAreaIdx);
  binningData.add(binGlobalIdx);

  return binningData;
}

atlas::FieldSet createGlobalHorizontalAverageWeights(const eckit::mpi::Comm & localComm,
                                                     const std::string & binType,
                                                     const std::size_t sizeOwned,
                                                     const atlas::FunctionSpace & fs) {
  atlas::Field gridAreaWgt =
    atlas::Field(binType + " weights",  atlas::array::make_datatype<double>(),
                 atlas::array::make_shape(1, sizeOwned));
  const std::string binTypeIdx = binType + " local indices";
  atlas::Field gridAreaIdx =
    atlas::Field(binTypeIdx, atlas::array::make_datatype<std::int32_t>(),
                 atlas::array::make_shape(1, sizeOwned));
  const std::string binTypeGIdx = binType + " global bins";
  atlas::Field binGlobalIdx =
    atlas::Field(binTypeGIdx, atlas::array::make_datatype<std::int32_t>(),
                 atlas::array::make_shape(1));

  atlas::Field weights;
  if (fs.type() == "StructuredColumns") {
    auto scFspace = atlas::functionspace::StructuredColumns(fs);
    std::string s = scFspace.grid().name().substr(0, 1);
    // "F" denotes a regular Gaussian grid.
    if (s == "F") weights = GaussCellArea(localComm, scFspace, sizeOwned, binType);
  } else if (fs.type() == "NodeColumns") {
    auto nodeFspace = atlas::functionspace::NodeColumns(fs);
    std::string s = nodeFspace.mesh().grid().name();
    if (s.substr(0, 6) == "CS-LFR") {
      weights = CSCellArea(localComm, nodeFspace, sizeOwned, binType);
    }
  }
  auto gridAreaWgtView = atlas::array::make_view<double, 2>(gridAreaWgt);
  auto gridAreaIdxView = atlas::array::make_view<std::int32_t, 2>(gridAreaIdx);
  auto weightsView = atlas::array::make_view<double, 2>(weights);
  for (std::size_t jn = 0; jn < sizeOwned; jn++) {
    gridAreaWgtView(0, jn) = weightsView(0, jn);
    gridAreaIdxView(0, jn) = static_cast<std::int32_t>(jn);
  }
  auto binGlobalIdxView = atlas::array::make_view<std::int32_t, 1>(binGlobalIdx);
  binGlobalIdxView(0) = 0;

  atlas::FieldSet binningData;
  binningData.add(gridAreaWgt);
  binningData.add(gridAreaIdx);
  binningData.add(binGlobalIdx);

  return binningData;
}

atlas::FieldSet createBinningData(const std::string & binType,
                                  const std::size_t sizeOwned,
                                  const oops::GeometryData & geomData) {
  atlas::FieldSet binningData;
  atlas::FunctionSpace fs = geomData.functionSpace();

  // populate area weights.
  atlas::FieldSet weightFlds;
  if (binType.compare("horizontal global average") == 0) {
    weightFlds = createGlobalHorizontalAverageWeights(geomData.comm(), binType,
                                                      sizeOwned, fs);
  } else if (binType.compare("horizontal grid point") == 0) {
    weightFlds = createHorizontalGridPointWeights(binType, sizeOwned, fs);
  } else {
    throw eckit::UserError(binType + " not accounted for when setting Binning Data", Here());
  }

  for (const atlas::Field & w : weightFlds) {
    binningData.add(w);
  }

  const atlas::Field lonlat = fs.lonlat();
  auto areaIdxView =
    atlas::array::make_view<const std::int32_t, 2>(binningData[binType + " local indices"]);
  auto lon = atlas::Field("longitude", atlas::array::make_datatype<double>(),
    atlas::array::make_shape(areaIdxView.shape()[0], areaIdxView.shape()[1]));
  auto lat = atlas::Field("latitude", atlas::array::make_datatype<double>(),
    atlas::array::make_shape(areaIdxView.shape()[0], areaIdxView.shape()[1]));
  auto lonlatView = atlas::array::make_view<const double, 2>(lonlat);
  auto latView = atlas::array::make_view<double, 2>(lat);
  auto lonView = atlas::array::make_view<double, 2>(lon);
  for (atlas::idx_t b = 0; b < areaIdxView.shape()[0]; ++b) {
    for (atlas::idx_t i = 0; i < areaIdxView.shape()[1]; ++i) {
      lonView(b, i) = lonlatView(areaIdxView(b, i), atlas::LON);
      latView(b, i) = lonlatView(areaIdxView(b, i), atlas::LAT);
    }
  }
  binningData.add(lat);
  binningData.add(lon);

  return binningData;
}

std::size_t getTotalBins(const eckit::mpi::Comm & localComm,
                         const std::string & binType,
                         const std::size_t sizeOwned) {
  std::size_t totalBins(0);
  if (binType.compare("horizontal global average") == 0) {
    totalBins = 1;
  } else if (binType.compare("horizontal grid point") == 0) {
    localComm.allReduce(sizeOwned, totalBins, eckit::mpi::sum());
  } else {
    throw eckit::UserError(binType + " not accounted for it getting TotalBins", Here());
  }
  return totalBins;
}

std::vector<std::size_t> getTotalPtsPerGlobalBin(const eckit::mpi::Comm & localComm,
                                                 const std::string & binType,
                                                 const atlas::FieldSet & binningData,
                                                 const std::size_t & totalBins) {
  std::vector<std::size_t> totalPtsPerGlobalBinLocal(totalBins, 0);
  std::vector<std::size_t> totalPtsPerGlobalBin(totalBins, 0);

  // access local indices
  std::string binTypeIdx = binType + " local indices";
  std::string binTypeGIdx = binType + " global bins";
  auto gridAreaIdxView = atlas::array::make_view<std::int32_t, 2>(binningData[binTypeIdx]);
  auto binGlobalIdxView = atlas::array::make_view<std::int32_t, 1>(binningData[binTypeGIdx]);

  for (atlas::idx_t b = 0; b < gridAreaIdxView.shape()[0]; b++) {
    // TODO(Marek) Note that in the general case this will not work
    // ... where the number of indices on each PE for each valid bin is not the same.
    // In that case we either extend code to use zero-padding and an addtional
    // binning field that states the maximum size for each bin OR
    // we extend the number of binning fields so that there is a separate array of each local bin
    // on each PE.
    totalPtsPerGlobalBinLocal[binGlobalIdxView(b)] = gridAreaIdxView.shape()[1];
  }

  localComm.allReduce(totalPtsPerGlobalBinLocal, totalPtsPerGlobalBin, eckit::mpi::sum());

  return totalPtsPerGlobalBin;
}

// If the calibration is switched on
// it will allocate fields (populated to zero) ready for covariance
// accumulation.
// TODO(Marek): maybe remove params here and just use netCDFConf_?
atlas::FieldSet createEnsembleStatsFSet(const std::string & binType,
                                        const WriteVariancesParameters & params,
                                        const eckit::LocalConfiguration & netCDFConf,
                                        const atlas::FieldSet & binningData) {
  atlas::FieldSet ensembleStats;
  if (params.calibrationParams.value() != boost::none) {
    oops::Variables fvars(params.fieldNames);
    std::size_t noOfStats = fvars.size();
    eckit::LocalConfiguration conf = params.toConfiguration();
    if (conf.has("additional cross covariances")) {
       noOfStats += conf.getSubConfigurations("additional cross covariances").size();
    }
    const std::string statsType = params.statisticsType;
    for (std::size_t s = 0; s < noOfStats; ++s) {
      atlas::Field fld;
      const std::string name = statsType + " " + std::to_string(s);
      const std::int32_t levels1 =
        util::getAttributeValue<std::int32_t>(netCDFConf, name, "levels 1");
      const std::int32_t levels2 =
        util::getAttributeValue<std::int32_t>(netCDFConf, name, "levels 2");
      const std::size_t bins = binningData[binType + " weights"].shape(0);
      if (statsType == "variance") {
        ASSERT(levels1 == levels2);
        fld = atlas::Field(name, atlas::array::make_datatype<double>(),
                           atlas::array::make_shape(bins, levels1));
        atlas::array::make_view<double, 2>(fld).assign(0.0);
        ensembleStats.add(fld);
      } else if (statsType == "vertical covariance") {
        fld = atlas::Field(name, atlas::array::make_datatype<double>(),
                           atlas::array::make_shape(bins, levels1, levels2));
        atlas::array::make_view<double, 3>(fld).assign(0.0);
        ensembleStats.add(fld);
      } else {
        throw eckit::UserError(statsType + " not accounted for", Here());
      }
    }
  }

  return ensembleStats;
}

}  // namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<WriteVariances>
  makerWriteVariances_("write variances");

// -----------------------------------------------------------------------------
void WriteVariances::writeInstantVariances(const eckit::mpi::Comm & comm,
                                           const atlas::FieldSet & fset,
                                           const std::string & description) const {
  // Select the fields to write out.
  std::size_t root(0);
  // If the `fieldNames` list is empty, write out all fields in `fset`.
  atlas::FieldSet fsetWrite;
  oops::Variables variablesToWrite;
  if (params_.fieldNames.value().empty()) {
    variablesToWrite = oops::Variables{fset.field_names()};
  } else {
    variablesToWrite = oops::Variables(params_.fieldNames.value());
  }
  for (const auto & field : fset) {
    if (variablesToWrite.has(field.name())) {
      fsetWrite.add(field);
    }
  }

  // Check there is at least one field to write out.
  if (fsetWrite.size() == 0) {
    oops::Log::warning() << "No field information to output." << std::endl;
    return;
  }

  const varianceInstantaneousParameters instParams =
    params_.instantaneousParams.value().value();
  std::stringstream filepath;
  filepath << instParams.outputPath.value() << "/"
           << description;
  const std::string filepathnc = filepath.str() + ".nc";

  if (params_.saveNetCDFFile) {
    eckit::LocalConfiguration conf;
    conf.set("filepath", filepath.str());

    // todo: replace the use of `stat` with `std::filesystem::exists`
    // when compilers have been upgraded.
    struct stat buffer;
    if (stat(filepathnc.c_str(), &buffer) == 0) {
      oops::Log::warning() << "File " << filepathnc
                           << " already exists. Overwriting." << std::endl;
    }

    // Write fsetWrite to file.
    const std::vector<std::string> dim_names{"bin index",
                                             "model levels"};
    const std::vector<atlas::idx_t> dim_sizes{fsetWrite[0].shape()[0],
                                              fsetWrite[0].shape()[1]};
    oops::Variables fsetVars(variablesToWrite);

    std::vector<std::vector<std::string>> dim_names_for_every_var;

    eckit::LocalConfiguration netcdfMetaData;
    for (const oops::Variable & var : fsetVars) {
      dim_names_for_every_var.push_back(dim_names);
      util::setAttribute<std::string>(
        netcdfMetaData, var.name(), "statistics type", "string",
        "horizontally-averaged variance");
    }

    std::vector<int> netcdf_general_ids;
    std::vector<int> netcdf_dim_ids;
    std::vector<int> netcdf_var_ids;
    std::vector<std::vector<int>> netcdf_dim_var_ids;

    if (comm.rank() == root) {
      ::util::atlasArrayWriteHeader(filepathnc,
                                    dim_names,
                                    dim_sizes,
                                    fsetVars,
                                    dim_names_for_every_var,
                                    netcdfMetaData,
                                    netcdf_general_ids,
                                    netcdf_dim_ids,
                                    netcdf_var_ids,
                                    netcdf_dim_var_ids);

      std::size_t t(0);
      for (const atlas::Field & fld : fsetWrite) {
        auto fview = atlas::array::make_view<const double, 2>(fld);
        ::util::atlasArrayWriteData(netcdf_general_ids,
                                    netcdf_var_ids[t],
                                    fview);
         ++t;
      }
      nc_close(netcdf_general_ids[0]);
    }

    // Output filename to test stream.
    oops::Log::test() << "Wrote file " << filepathnc << std::endl;
  } else {
    // Output filename to test stream.
    oops::Log::test() << "Did not write file " << filepathnc << std::endl;
  }
}

// currently only does "horizontally-averaged variances"
void WriteVariances::diagnostics(const std::string & tag,
                                 const oops::FieldSet3D & fset) const {
  atlas::FieldSet varianceFset;
  const std::size_t nbins = 1;
  const std::size_t root = 0;
  computeVarianceFieldSetInstant(innerGeometryData_.comm(), tag,
                                 fset.fieldSet(), root, nbins,
                                 varianceFset);
  // note that if this is to be used in a parallel farming context.
  // we will need both the global root, local root and both global and local
  // eckit Communicators.
  printInstantBinnedVariances(innerGeometryData_.comm(), root,
                              tag, nbins, varianceFset);
  writeInstantVariances(innerGeometryData_.comm(), varianceFset, tag);
}

// -----------------------------------------------------------------------------

WriteVariances::WriteVariances(const oops::GeometryData & outerGeometryData,
                               const oops::Variables & outerVars,
                               const eckit::Configuration & covarConfig,
                               const Parameters_ & params,
                               const oops::FieldSet3D & xb,
                               const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    innerVars_(outerVars),
    params_(params),
    binType_(getBinType(params_)),
    netCDFConf_(createNetCDFHeaderInput(binType_, params, outerVars)),
    sizeOwned_(getSizeOwned(innerGeometryData_.functionSpace())),
    binningData_(createBinningData(binType_, sizeOwned_, innerGeometryData_)),
    totalBins_(getTotalBins(innerGeometryData_.comm(), binType_, sizeOwned_)),
    totalPtsPerGlobalBin_(
      getTotalPtsPerGlobalBin(innerGeometryData_.comm(), binType_,
                              binningData_, totalBins_)),
    ensembleStats_(
      createEnsembleStatsFSet(binType_, params_, netCDFConf_, binningData_)),
    datetime_(xb.validTime()),
    count_multiply_(1),
    count_multiplyad_(1),
    count_leftinversemultiply_(1),
    sample_size_(0)
{
  oops::Log::trace() << classname() << "::WriteVariances starting" << std::endl;

  const binningParameters bparams = params.binning;
  // By setting binning filePath and binning mpiRankPattern we can
  // optionally dump the binning fields to file of each mpi rank.
  if ((bparams.filePath.value() != ::boost::none) &&
      (bparams.mpiRankPattern.value() != ::boost::none)) {
    const std::string mpi_pattern = bparams.mpiRankPattern.value().value();
    const std::string mpi_rank = std::to_string(eckit::mpi::comm().rank());
    eckit::LocalConfiguration binConfig;
    bparams.serialize(binConfig);
    ::util::seekAndReplace(binConfig, mpi_pattern, mpi_rank);
    const std::string ncfilepath = "./" + binConfig.getString("file path");

    // create a processed binningData local to each PE.
    atlas::FieldSet processedBinnedData;

    // TODO(Marek) extend to addtional bin types and dump a bin file for each rank
    if (((binType_ == "horizontal global average") ||
        (binType_ == "horizontal grid point")) &&
         bparams.oneFilePerTask) {
      processedBinnedData.add(binningData_["longitude"]);
      processedBinnedData.add(binningData_["latitude"]);
      processedBinnedData.add(binningData_[binType_ + " weights"]);
      processedBinnedData.add(binningData_[binType_ + " global bins"]);


      // write to file
      const std::string gBinTypeIdx = "global PE bin index";
      const std::string binTypeIdx = "local PE bin index";
      const std::string levelsIdx1 = "local PE horizontal point index";

      const std::vector<std::string> dim_names{gBinTypeIdx, binTypeIdx, levelsIdx1};
      const std::vector<atlas::idx_t> dim_sizes{processedBinnedData[0].shape()[0],
                                                processedBinnedData[0].shape()[0],
                                                processedBinnedData[0].shape()[1]};

      // adding global header
      eckit::LocalConfiguration netCDFConf = netCDFConf_;
      util::setAttribute<std::string>(
        netCDFConf, "global metadata", "date time", "string", datetime_.toString());
      util::setAttribute<std::string>(
        netCDFConf, "global metadata", "binning type", "string", binType_);

      oops::Variables fsetVars;
      for (const std::string & name : processedBinnedData.field_names()) {
        oops::VariableMetaData varMeta =
          name == binType_ + " global bins" ?
          oops::VariableMetaData(oops::VerticalStagger::CENTER, oops::ModelDataType::Int32) :
          oops::VariableMetaData(oops::VerticalStagger::CENTER, oops::ModelDataType::Real64);
        fsetVars.push_back(oops::Variable(name, varMeta));
      }

      std::vector<int> netcdf_general_ids;
      std::vector<int> netcdf_dim_ids;
      std::vector<int> netcdf_var_ids;
      std::vector<std::vector<int>> netcdf_dim_var_ids;
      std::vector<std::vector<std::string>>
          dim_names_for_every_var;
      for (auto & f : fsetVars) {
        if (f.name() == binType_ + " global bins") {
          dim_names_for_every_var.push_back(
            std::vector<std::string>{dim_names[0]});
        } else {
          dim_names_for_every_var.push_back(
            std::vector<std::string>{dim_names[1], dim_names[2]});
        }
      }

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

      for (const atlas::Field & fld : processedBinnedData) {
        if (fld.shape().size() == 2) {
          auto fview = atlas::array::make_view<const double, 2>(fld);
          ::util::atlasArrayWriteData(netcdf_general_ids,
                                      netcdf_var_ids[t],
                                      fview);
        } else if (fld.shape().size() == 1) {
            auto fview = atlas::array::make_view<const std::int32_t, 1>(fld);
            ::util::atlasArrayWriteData(netcdf_general_ids,
                                        netcdf_var_ids[t],
                                        fview);
        }
        ++t;
      }

      int retval;
      if ((retval = nc_close(netcdf_general_ids[0]))) {
        throw eckit::Exception("NetCDF closing error", Here());
      }
    }
  }
  oops::Log::trace() << classname() << "::WriteVariances done" << std::endl;
}

// -----------------------------------------------------------------------------

void WriteVariances::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  if (params_.instantaneousParams.value() != ::boost::none) {
    const varianceInstantaneousParameters instParams =
      params_.instantaneousParams.value().value();
    if (instParams.multiplyFileName.value() != ::boost::none) {
      const std::string tag = instParams.multiplyFileName.value().value() +
        "_" + std::to_string(count_multiply_++);
      diagnostics(tag, fset);
    }
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void WriteVariances::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  if (params_.instantaneousParams.value() != ::boost::none) {
    const varianceInstantaneousParameters instParams =
      params_.instantaneousParams.value().value();
    if (instParams.multiplyADFileName.value() != ::boost::none) {
      const std::string tag = instParams.multiplyADFileName.value().value() +
        "_" + std::to_string(count_multiplyad_++);
      diagnostics(tag, fset);
    }
  }

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void WriteVariances::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  if (params_.instantaneousParams.value() != ::boost::none) {
    const varianceInstantaneousParameters instParams =
      params_.instantaneousParams.value().value();
    if (instParams.leftInverseFileName.value() != ::boost::none) {
      const std::string tag = instParams.leftInverseFileName.value().value() +
        "_" + std::to_string(count_leftinversemultiply_++);
      diagnostics(tag, fset);
    }
  }

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void WriteVariances::directCalibration(const oops::FieldSets & fset) {
  // loop over ensemble members in fset to calculate statistics
  // (dependent on statistics type)
  const std::string statsType = params_.statisticsType;
  if (fset.ens_size() > 0) {
    if (statsType == "vertical covariance") {
      sample_size_ = updateVerticalCovariances(
        netCDFConf_, binningData_, fset, sample_size_, ensembleStats_);
    } else if (statsType == "variance") {
      sample_size_ = updateVariances(
        netCDFConf_, binningData_, fset, sample_size_, ensembleStats_);
    }
  }
}

// -----------------------------------------------------------------------------

void WriteVariances::write() const {
  const util::calibrationWriteParameters writeParams =
    (params_.calibrationParams.value())->writeParams.value();
  eckit::LocalConfiguration writeConfig;
  writeParams.serialize(writeConfig);

  const std::string mpi_pattern = writeParams.mpiPattern;
  const std::string mpi_size = std::to_string(eckit::mpi::comm().size());
  ::util::seekAndReplace(writeConfig, mpi_pattern, mpi_size);
  const std::string ncfilepath = "./" + writeConfig.getString("file path");

  using atlas::array::make_view;
  using atlas::idx_t;

  // gather and sum on pe 0
  const std::size_t root(0);
  const std::string binType =
    util::getAttributeValue<std::string>(netCDFConf_,
                                         ensembleStats_[0].name(),
                                         "binning type");
  atlas::FieldSet ensembleStatsToWrite =
    util::gatherSumFieldSet(innerGeometryData_.comm(),
                            root,
                            totalBins_,
                            binningData_[binType + " global bins"],
                            ensembleStats_);

  const std::string binTypeIdx = binType + " index";
  const std::string levelsIdx1 = "levels index 1";
  const std::string levelsIdx2 = "levels index 2";

  std::vector<std::string> dim_names{binTypeIdx, levelsIdx1};
  std::vector<idx_t> dim_sizes{ensembleStatsToWrite[0].shape()[0],
                               ensembleStatsToWrite[0].shape()[1]};
  if (ensembleStatsToWrite[0].shape().size() == 3) {
    dim_names.push_back(levelsIdx2);
    dim_sizes.push_back(ensembleStatsToWrite[0].shape()[1]);
  }

  // adding global header
  eckit::LocalConfiguration netCDFConf = netCDFConf_;
  util::setAttribute<std::string>(
    netCDFConf, "global metadata", "covariance name", "string", writeParams.covName);
  util::setAttribute<std::string>(
    netCDFConf, "global metadata", "date time", "string", datetime_.toString());
  util::setAttribute<std::int32_t>(
    netCDFConf, "global metadata", "no of samples", "int32",
    static_cast<std::int32_t>(sample_size_));

  // Note will need to change fsetVars to allow for lat, lon index and weight
  const oops::Variables fsetVars(ensembleStatsToWrite.field_names());

  std::vector<int> netcdf_general_ids;
  std::vector<int> netcdf_dim_ids;
  std::vector<int> netcdf_var_ids;
  std::vector<std::vector<int>> netcdf_dim_var_ids;
  std::vector<std::vector<std::string>>
    dim_names_for_every_var(fsetVars.size(), dim_names);

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

    if (ensembleStatsToWrite[0].shape().size() == 2) {
      for (const atlas::Field & fld : ensembleStatsToWrite) {
        auto fview = make_view<const double, 2>(fld);
        ::util::atlasArrayWriteData(netcdf_general_ids,
                                    netcdf_var_ids[t],
                                    fview);
        ++t;
      }
    } else if (ensembleStatsToWrite[0].shape().size() == 3) {
      for (const atlas::Field & fld : ensembleStatsToWrite) {
        auto fview = make_view<const double, 3>(fld);
        ::util::atlasArrayWriteData(netcdf_general_ids,
                                    netcdf_var_ids[t],
                                    fview);
         ++t;
      }
    }

    int retval;
    if ((retval = nc_close(netcdf_general_ids[0]))) {
      throw eckit::Exception("NetCDF closing error", Here());
    }
  }
}

// -----------------------------------------------------------------------------

void WriteVariances::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
