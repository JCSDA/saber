/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/WriteVariances.h"

#include <netcdf.h>
#include <sys/stat.h>

#include <algorithm>
#include <set>
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
    const atlas::Grid & grid = atlas::functionspace::NodeColumns(fs).mesh().grid();
    if (grid) {
       const std::string & gridname = grid.name();
       if (gridname.substr(0, 6) == "CS-LFR") {
         sizeOwned = atlas::functionspace::CubedSphereNodeColumns(fs).sizeOwned();
         fsViable = true;
       }
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
  util::setAttribute<std::string>(netCDFConf, binWeights, "binning_type",
                                  "string", binType);
  std::string binIndices = binType + " indices";
  util::setAttribute<std::string>(netCDFConf, binIndices, "long_name",
                                  "string", binIndices);
  util::setAttribute<std::string>(netCDFConf, binIndices, "binning_type",
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
  const std::string binTypeHExtent = binType + " horizontal extent";
  atlas::Field binHorizExtent =
    atlas::Field(binTypeHExtent, atlas::array::make_datatype<std::int32_t>(),
                 atlas::array::make_shape(sizeOwned));

  auto gridAreaWgtView = atlas::array::make_view<double, 2>(gridAreaWgt);
  auto gridAreaIdxView = atlas::array::make_view<std::int32_t, 2>(gridAreaIdx);
  auto globalIdxView = atlas::array::make_view<atlas::gidx_t, 1>(fs.global_index());
  auto binGlobalIdxView = atlas::array::make_view<std::int32_t, 1>(binGlobalIdx);
  auto binHorizExtentView = atlas::array::make_view<std::int32_t, 1>(binHorizExtent);

  for (std::size_t jn = 0; jn < sizeOwned; jn++) {
    gridAreaWgtView(jn, 0) = 1.0;
    gridAreaIdxView(jn, 0) = static_cast<std::int32_t>(jn);
    binGlobalIdxView(jn) = static_cast<std::int32_t>(globalIdxView(jn)-1);
    binHorizExtentView(jn) = 1.0;
  }
  atlas::FieldSet binningData;
  binningData.add(gridAreaWgt);
  binningData.add(gridAreaIdx);
  binningData.add(binGlobalIdx);
  binningData.add(binHorizExtent);

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
  const std::string binTypeHExtent = binType + " horizontal extent";
  atlas::Field binHorizExtent =
    atlas::Field(binTypeHExtent, atlas::array::make_datatype<std::int32_t>(),
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
  auto binHorizExtentView = atlas::array::make_view<std::int32_t, 1>(binHorizExtent);
  binHorizExtentView(0) = sizeOwned;

  atlas::FieldSet binningData;
  binningData.add(gridAreaWgt);
  binningData.add(gridAreaIdx);
  binningData.add(binGlobalIdx);
  binningData.add(binHorizExtent);

  return binningData;
}

// Idea here is to initially
atlas::FieldSet createOverlappingLatitudeWeights(const eckit::mpi::Comm & localComm,
                                                 const std::string & binType,
                                                 const std::size_t sizeOwned,
                                                 const std::size_t totalBins,
                                                 const atlas::FunctionSpace & fs) {
  atlas::FieldSet globalAverageWeights =
    createGlobalHorizontalAverageWeights(localComm,
                                         std::string{"horizontal global average"},
                                         sizeOwned, fs);
  atlas::idx_t asizeOwned = sizeOwned;
  std::vector<double> latMin(totalBins, 0.0);
  std::vector<double> latMax(totalBins, 0.0);
  std::vector<double> latPeak(totalBins, 0.0);
  double deltaLat = 180.0 / static_cast<double>(totalBins);
  for (std::size_t b = 0; b < totalBins; ++b) {
    latPeak[b] = -90.0 + (b + 0.5) * deltaLat;
    latMin[b] = -90.0 + (b - 0.5) * deltaLat;
    latMax[b] = -90.0 + (b + 1.5) * deltaLat;
  }
  latMin[0] = -90.0;
  latMin[1] = latPeak[0];
  latMax[totalBins-1] = 90.0;
  latMax[totalBins-2] = latPeak[totalBins-1];

  const atlas::Field lonlat = fs.lonlat();
  auto lonlatView = atlas::array::make_view<const double, 2>(lonlat);

  // need ideally to have no of local bins on each pe.
  // need no of local horizontal pts on each pe for each bin.
  std::set<std::size_t> globalBins;
  for (atlas::idx_t jn = 0; jn < asizeOwned; ++jn) {
    for (std::size_t b = 0; b < totalBins; ++b) {
      if (lonlatView(jn, atlas::LAT) >= latMin[b] && lonlatView(jn, atlas::LAT) <= latMax[b]) {
        globalBins.insert(b);
      }
    }
  }

  std::vector<std::size_t> horizontalIndexExtent(globalBins.size());

  // total weight should be for all the bins
  std::vector<double> totalweights(totalBins, 0.0);


  auto weightsGAverageView = atlas::array::make_view<double, 2>(
    globalAverageWeights[std::string{"horizontal global average weights"}]);

  std::size_t localBinIndx(0);
  for (std::size_t b : globalBins) {
    std::size_t extent(0);
    for (atlas::idx_t jn = 0; jn < asizeOwned; ++jn) {
      if (lonlatView(jn, atlas::LAT) >= latMin[b] && lonlatView(jn, atlas::LAT) <= latMax[b]) {
        double latweight(0.0);
        if (((b == 0) && (lonlatView(jn, atlas::LAT) <= latPeak[b])) ||
            ((b == totalBins - 1) && (lonlatView(jn, atlas::LAT) >= latPeak[b]))) {
          latweight = 1.0;
        } else {
          latweight = 1.0 - std::abs(lonlatView(jn, atlas::LAT) - latPeak[b]) / deltaLat;
        }
        ASSERT(latweight >= 0);
        totalweights[b] += weightsGAverageView(0, jn) * latweight;
        ++extent;
      }
    }
    horizontalIndexExtent[localBinIndx] = extent;
    localBinIndx++;
  }
  // need a sum totalweights across PEs
  localComm.allReduceInPlace(totalweights.begin(), totalweights.end(), eckit::mpi::sum());

  auto maxElement = std::max_element(horizontalIndexExtent.begin(), horizontalIndexExtent.end());
  std::size_t maxExtent = *maxElement;

  atlas::Field gridAreaWgt =
    atlas::Field(binType + " weights",  atlas::array::make_datatype<double>(),
                 atlas::array::make_shape(globalBins.size(), maxExtent));
  const std::string binTypeIdx = binType + " local indices";
  atlas::Field gridAreaIdx =
    atlas::Field(binTypeIdx, atlas::array::make_datatype<std::int32_t>(),
                 atlas::array::make_shape(globalBins.size(), maxExtent));
  const std::string binTypeGIdx = binType + " global bins";
  atlas::Field binGlobalIdx =
    atlas::Field(binTypeGIdx, atlas::array::make_datatype<std::int32_t>(),
                 atlas::array::make_shape(globalBins.size()));
  const std::string binTypeHExtent = binType + " horizontal extent";
  atlas::Field binHorizExtent =
    atlas::Field(binTypeHExtent, atlas::array::make_datatype<std::int32_t>(),
                 atlas::array::make_shape(globalBins.size()));

  auto gridAreaWgtView = atlas::array::make_view<double, 2>(gridAreaWgt);
  auto gridAreaIdxView = atlas::array::make_view<std::int32_t, 2>(gridAreaIdx);
  auto binGlobalIdxView = atlas::array::make_view<std::int32_t, 1>(binGlobalIdx);
  auto binHorizExtentView = atlas::array::make_view<std::int32_t, 1>(binHorizExtent);
  gridAreaWgtView.assign(0.0);
  gridAreaIdxView.assign(0.0);
  binGlobalIdxView.assign(0.0);
  binHorizExtentView.assign(0.0);

  localBinIndx = 0;
  for (std::size_t b : globalBins) {
    binGlobalIdxView(localBinIndx) = b;
    binHorizExtentView(localBinIndx) = horizontalIndexExtent[localBinIndx];

    std::size_t indx(0);
    for (atlas::idx_t jn = 0; jn < asizeOwned; ++jn) {
      if (lonlatView(jn, atlas::LAT) >= latMin[b] && lonlatView(jn, atlas::LAT) <= latMax[b]) {
        double latweight(0.0);
        if (((b == 0) && (lonlatView(jn, atlas::LAT) <= latPeak[b])) ||
            ((b == totalBins - 1) && (lonlatView(jn, atlas::LAT) >= latPeak[b]))) {
          latweight = 1.0;
        } else {
          latweight = 1.0 - std::abs(lonlatView(jn, atlas::LAT) - latPeak[b]) / deltaLat;
        }

        gridAreaIdxView(localBinIndx, indx) = jn;
        gridAreaWgtView(localBinIndx, indx) =
          weightsGAverageView(0, jn) * latweight / totalweights[b];
        ++indx;
      }
    }
    localBinIndx++;
  }
  atlas::FieldSet binningData;
  binningData.add(gridAreaWgt);
  binningData.add(gridAreaIdx);
  binningData.add(binGlobalIdx);
  binningData.add(binHorizExtent);

  return binningData;
}

atlas::FieldSet createBinningData(const std::string & binType,
                                  const std::size_t sizeOwned,
                                  const std::size_t noOfBins,
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
  } else if (binType.compare("overlapping area-weighted latitude bands") == 0) {
    weightFlds = createOverlappingLatitudeWeights(geomData.comm(), binType,
                                                  sizeOwned, noOfBins, fs);
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
                         const std::size_t sizeOwned,
                         const binningParameters & bparams) {
  std::size_t totalBins(0);
  if (binType.compare("horizontal global average") == 0) {
    totalBins = 1;
  } else if (binType.compare("horizontal grid point") == 0) {
    localComm.allReduce(sizeOwned, totalBins, eckit::mpi::sum());
  } else if (binType.compare("overlapping area-weighted latitude bands") == 0) {
    if (bparams.noOfBins.value() == ::boost::none) {
      throw eckit::UserError("no of bins needs to be set for when using " + binType , Here());
    }
    totalBins = bparams.noOfBins.value().value();
  } else {
    throw eckit::UserError(binType + " not accounted for it getting TotalBins", Here());
  }
  return totalBins;
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
      std::string netCDFShortName = statsType;
      std::replace(netCDFShortName.begin(), netCDFShortName.end(), ' ', '_');
      netCDFShortName.append("_" + std::to_string(s));

      const std::int32_t levels1 =
        util::getAttributeValue<std::int32_t>(netCDFConf, netCDFShortName, "levels_1");
      const std::int32_t levels2 =
        util::getAttributeValue<std::int32_t>(netCDFConf, netCDFShortName, "levels_2");
      const std::size_t bins = binningData[binType + " weights"].shape(0);
      if (statsType == "variance") {
        ASSERT(levels1 == levels2);
        fld = atlas::Field(netCDFShortName, atlas::array::make_datatype<double>(),
                           atlas::array::make_shape(bins, levels1));
        atlas::array::make_view<double, 2>(fld).assign(0.0);
        ensembleStats.add(fld);
      } else if (statsType == "vertical covariance") {
        fld = atlas::Field(netCDFShortName, atlas::array::make_datatype<double>(),
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
    const std::vector<std::string> dim_names{"bin_index",
                                             "model_levels"};
    const std::vector<atlas::idx_t> dim_sizes{fsetWrite[0].shape()[0],
                                              fsetWrite[0].shape()[1]};
    oops::Variables fsetVars(variablesToWrite);

    std::vector<std::vector<std::string>> dim_names_for_every_var;

    eckit::LocalConfiguration netcdfMetaData;
    for (const oops::Variable & var : fsetVars) {
      dim_names_for_every_var.push_back(dim_names);
      util::setAttribute<std::string>(
        netcdfMetaData, var.name(), "statistics_type", "string",
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
    totalBins_(getTotalBins(innerGeometryData_.comm(), binType_, sizeOwned_, params_.binning)),
    binningData_(createBinningData(binType_, sizeOwned_, totalBins_, innerGeometryData_)),
    ensembleStats_(
      createEnsembleStatsFSet(binType_, params_, netCDFConf_, binningData_)),
    datetime_(xb.validTime()),
    count_multiply_(1),
    count_multiplyad_(1),
    count_leftinversemultiply_(1),
    sample_size_(0)
{
  oops::Log::info() << classname() << "::WriteVariances starting" << std::endl;

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
        (binType_ == "horizontal grid point") ||
        (binType_ == "overlapping area-weighted latitude bands")) &&
         bparams.oneFilePerTask) {
      processedBinnedData.add(binningData_["longitude"]);
      processedBinnedData.add(binningData_["latitude"]);
      processedBinnedData.add(binningData_[binType_ + " weights"]);
      processedBinnedData.add(binningData_[binType_ + " global bins"]);
      processedBinnedData.add(binningData_[binType_ + " horizontal extent"]);

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
        netCDFConf, "global metadata", "date_time", "string", datetime_.toString());
      util::setAttribute<std::string>(
        netCDFConf, "global metadata", "binning_type", "string", binType_);

      oops::Variables fsetVars;
      for (const std::string & name : processedBinnedData.field_names()) {
        oops::VariableMetaData varMeta;
        if ( name == binType_ + " global bins" ) {
          varMeta =
            oops::VariableMetaData(oops::VerticalStagger::CENTER, oops::ModelDataType::Int32);
        } else {
          varMeta =
            oops::VariableMetaData(oops::VerticalStagger::CENTER, oops::ModelDataType::Real64);
        }
        fsetVars.push_back(oops::Variable(name, varMeta));
      }

      std::vector<int> netcdf_general_ids;
      std::vector<int> netcdf_dim_ids;
      std::vector<int> netcdf_var_ids;
      std::vector<std::vector<int>> netcdf_dim_var_ids;
      std::vector<std::vector<std::string>>
          dim_names_for_every_var;
      for (auto & f : fsetVars) {
        if ((f.name() == binType_ + " global bins") ||
            (f.name() == binType_ + " horizontal extent"))
        {
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
  oops::Log::info() << classname() << "::WriteVariances done" << std::endl;
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
                                         "binning_type");
  atlas::FieldSet ensembleStatsToWrite =
    util::gatherSumFieldSet(innerGeometryData_.comm(),
                            root,
                            totalBins_,
                            binningData_[binType + " global bins"],
                            ensembleStats_);

  const std::string binTypeIdx = "binning_index";
  const std::string levelsIdx1 = "levels_index_1";
  const std::string levelsIdx2 = "levels_index_2";

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
    netCDFConf, "global metadata", "covariance_name", "string", writeParams.covName);
  util::setAttribute<std::string>(
    netCDFConf, "global metadata", "date_time", "string", datetime_.toString());
  util::setAttribute<std::int32_t>(
    netCDFConf, "global metadata", "no_of_samples", "int32",
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
