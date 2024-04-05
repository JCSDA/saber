/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/WriteVariances.h"

#include <sys/stat.h>

#include <sstream>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/util/Earth.h"
#include "atlas/util/GaussianLatitudes.h"

#include "eckit/mpi/Comm.h"

#include "oops/base/FieldSet3D.h"
#include "oops/util/AtlasArrayUtil.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace saber {
namespace vader {

namespace {

atlas::Field CSCellArea(const atlas::FunctionSpace& functionspace,
                        const std::string& fieldName)
{
  // Note that this does not work for CS cells on the North and South Poles.
  constexpr auto degreesToRadians = M_PI / 180.;
  const auto nodeColumns = atlas::functionspace::NodeColumns(functionspace);
  const auto cubedSphereGrid = atlas::CubedSphereGrid(nodeColumns.mesh().grid());
  const auto lonlat = atlas::array::make_view<double, 2>(functionspace.lonlat());
  const auto gridRes = cubedSphereGrid.N();
  constexpr double tolerance = 1e-8;

  const auto& proj = cubedSphereGrid.cubedSphereProjection();

  const double scaledCSAlphaBetaTileCellArea = M_PI/(2*gridRes) * M_PI/(2*gridRes) *
    atlas::util::Earth::radius() * atlas::util::Earth::radius();

  auto gridAreaField = functionspace.createField<double>(atlas::option::name(fieldName) |
                                                         atlas::option::levels(1));
  auto gridAreaView = atlas::array::make_view<double, 2>(gridAreaField);

  // Area is calculated by transforming the grid area of a cell on a cube onto the
  // Mercator projection and then converting it onto the sphere by applying the cosine.
  // This is only an approximation to the cell area and only converges linearly as a function
  // of resolution.
  for (std::size_t i = 0; i < gridAreaView.size(); i++) {
    const atlas::PointLonLat ll = atlas::PointLonLat(lonlat(i, atlas::LON),
                                                     lonlat(i, atlas::LAT));
    if (std::abs(std::abs(ll.lat()) - 90.0) < tolerance) {
      gridAreaView(i, 0) = scaledCSAlphaBetaTileCellArea;
    } else {
      const double cosLat = std::cos(degreesToRadians * ll.lat());
      const double gridJacDet = 1.0/proj.jacobian(ll).determinant();
      gridAreaView(i, 0) = cosLat * gridJacDet * scaledCSAlphaBetaTileCellArea;
    }
  }

  return gridAreaField;
}

atlas::Field GaussCellArea(atlas::functionspace::StructuredColumns & structCols,
                           const std::string& fieldName)
{
  const auto structGrid = atlas::StructuredGrid(structCols.grid());
  const double surfaceAreaOfEarth = 4.0 * M_PI * atlas::util::Earth::radius() *
    atlas::util::Earth::radius();
  auto lonlatView = atlas::array::make_view<double, 2>(structCols.lonlat());
  auto northSouthPts = structGrid.ny();
  std::vector<double> latitudes(northSouthPts, 0.0);
  std::vector<double> weights(northSouthPts, 0.0);
  // The function below calculates the Gaussian latitudes (latitudes) and a
  // normalized surface area assigned to each Gaussian latitude (weights)
  // The normalized surface area weights sum to 1.
  atlas::util::gaussian_quadrature_npole_spole(northSouthPts/2,
                                               latitudes.data(),
                                               weights.data());
  atlas::Field gridAreaField =
   structCols.createField<double>(
     atlas::option::name(fieldName) | atlas::option::levels(1));
  auto gridAreaView = atlas::array::make_view<double, 2>(gridAreaField);
  gridAreaView.assign(0.0);

  for (std::size_t i = 0; i < gridAreaView.size(); i++) {
    auto ll = atlas::PointLonLat(lonlatView(i, atlas::LON), lonlatView(i, atlas::LAT));

    // find index of lat
    for (std::size_t j = 0;  j < latitudes.size(); j++) {
      if (std::abs(latitudes[j] - ll.lat()) < __FLT_EPSILON__) {
        gridAreaView(i, 0) = weights[j] * surfaceAreaOfEarth /
          static_cast<double>(structGrid.nx(j));
        break;
      }
    }
  }
  return gridAreaField;
}

void computeBinnedVariances(const eckit::mpi::Comm & localComm,
                            const std::string & tag,
                            const atlas::Field & gridptAreaWeights,
                            const atlas::idx_t & sizeOwned,
                            const atlas::FieldSet & fset,
                            const std::size_t & localRoot,
                            std::size_t & nbins,
                            atlas::FieldSet & writeFset)
{
  nbins = 1;  // as we have a global variance for each field on each model level
  auto gridAreaView = atlas::array::make_view<double, 2>(gridptAreaWeights);
  double surfaceAreaOfEarth(0.0);
  for (atlas::idx_t jn = 0; jn < sizeOwned; ++jn) {
    surfaceAreaOfEarth += gridAreaView(jn, 0);
  }

  localComm.allReduceInPlace(surfaceAreaOfEarth, eckit::mpi::sum());

  for (const auto & field : fset) {
    auto fldView = atlas::array::make_view<double, 2>(fset[field.name()]);
    auto fldVariance = atlas::Field(field.name(), atlas::array::make_datatype<double>(),
      atlas::array::make_shape(nbins, fldView.shape(1)));
    auto fldVarianceView = atlas::array::make_view<double, 2>(fldVariance);
    fldVarianceView.assign(0.0);

    for (atlas::idx_t jl = 0; jl < fldView.shape(1); ++jl) {
      for (atlas::idx_t jn = 0; jn < sizeOwned; ++jn) {
        fldVarianceView(0, jl) += fldView(jn, jl) * fldView(jn, jl) * gridAreaView(jn, 0);
      }
      fldVarianceView(0, jl) /= surfaceAreaOfEarth;
    }

    // Note that the writeFset will only have the correct value on the root PE.
    util::gatherSum(localComm, localRoot, fldVarianceView);

    writeFset.add(fldVariance);
  }
}

// This should work for both cubed-sphere and gauss meshes
void printBinnedVariances(const eckit::mpi::Comm & globalComm,
                          const std::size_t & globalRoot,
                          const std::string & tag,
                          const std::size_t & nbins,
                          const atlas::FieldSet & writeFset)
{
  // print variance for each field
  for (const auto & f : writeFset) {
    auto fldVarianceView = atlas::array::make_view<double, 2>(f);
    std::vector<double> variance(nbins, 0.0);
    for (std::size_t b = 0; b < nbins; ++b) {
      for (atlas::idx_t l = 0; l < fldVarianceView.shape(1); ++l) {
        variance[b] += fldVarianceView(b, l);
      }
      variance[b] /= static_cast<double>(fldVarianceView.shape(1));
    }

    for (std::size_t b = 0; b < nbins; ++b) {
      oops::Log::test() << tag
                        << " ; name = " << f.name()
                        << " ; bin index = " << b
                        << " ; binned variance = " << variance[b]
                        << std::endl;
    }
  }
}

void computeVarianceFieldSet(const eckit::mpi::Comm & comm,
                             const std::string & tag,
                             const atlas::FieldSet & fset,
                             const std::size_t & root,
                             std::size_t & nbins,
                             atlas::FieldSet & varianceFset) {
  nbins = 1;
  bool computeVariance(false);
  if (fset[0].functionspace().type() == "StructuredColumns") {
    auto scFspace = atlas::functionspace::StructuredColumns(fset[0].functionspace());
    std::string s = scFspace.grid().name().substr(0, 1);
    // "F" denotes a regular Gaussian grid.
    if (s == "F") {
      atlas::Field gridptAreaWeights = GaussCellArea(scFspace,
                                                     "grid point area weights");
      computeBinnedVariances(comm,
                             tag,
                             gridptAreaWeights,
                             static_cast<atlas::idx_t>(scFspace.sizeOwned()),
                             fset,
                             root,
                             nbins,
                             varianceFset);
      computeVariance = true;
    }
  } else if (fset[0].functionspace().type() == "NodeColumns") {
    std::string s = atlas::functionspace::NodeColumns(
      fset[0].functionspace()).mesh().grid().name();
    if (s.substr(0, 6) == "CS-LFR") {
      auto csFspace = atlas::functionspace::CubedSphereNodeColumns(fset[0].functionspace());
      atlas::Field gridptAreaWeights = CSCellArea(csFspace,
                                                  "grid point area weights");
      computeBinnedVariances(comm,
                             tag,
                             gridptAreaWeights,
                             static_cast<atlas::idx_t>(csFspace.sizeOwned()),
                             fset,
                             root,
                             nbins,
                             varianceFset);
      computeVariance = true;
    }
  }
  if (!computeVariance) {
    oops::Log::error() << "WriteVariances:: variances calculations not supported "
                       << "for this specific mesh" << std::endl;
  }
}
}  // namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<WriteVariances>
  makerWriteVariances_("write variances");

// -----------------------------------------------------------------------------
void WriteVariances::writeToFile(const eckit::mpi::Comm & comm,
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

  // Produce a filename and increment the relevant counter.
  std::stringstream filepath;
  filepath << params_.outputPath.value() << "/"
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
    std::vector<std::string> field_names = fsetWrite.field_names();

    std::vector<std::vector<std::string>> dim_names_for_every_var;
    for (auto & field : field_names) {
      field.append(" horizontally-averaged variance");
      dim_names_for_every_var.push_back(dim_names);
    }

    std::vector<int> netcdf_general_ids;
    std::vector<int> netcdf_dim_ids;
    std::vector<int> netcdf_var_ids;
    std::vector<std::vector<int>> netcdf_dim_var_ids;

    if (comm.rank() == root) {
      ::util::atlasArrayWriteHeader(filepathnc,
                                    dim_names,
                                    dim_sizes,
                                    field_names,
                                    dim_names_for_every_var,
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

void WriteVariances::diagnostics(const std::string & tag,
                                 const oops::FieldSet3D & fset) const {
  atlas::FieldSet varianceFset;
  std::size_t nbins;
  computeVarianceFieldSet(innerGeometryData_.comm(), tag,
                          fset.fieldSet(), 0, nbins,
                          varianceFset);
  // note that if this is to be used in a parallel farming context.
  // we will need both the global root, local root and both global and local
  // eckit Communicators.
  printBinnedVariances(innerGeometryData_.comm(), 0,
                     tag, nbins, varianceFset);
  writeToFile(innerGeometryData_.comm(), varianceFset, tag);
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
    count_multiply_(1),
    count_multiplyad_(1),
    count_leftinversemultiply_(1)
{
  oops::Log::trace() << classname() << "::WriteVariances starting" << std::endl;

  oops::Log::trace() << classname() << "::WriteVariances done" << std::endl;
}

// -----------------------------------------------------------------------------

void WriteVariances::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  if (params_.multiplyFileName.value() != ::boost::none) {
    const std::string tag = params_.multiplyFileName.value().value() +
      "_" + std::to_string(count_multiply_++);
    diagnostics(tag, fset);
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void WriteVariances::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  if (params_.multiplyADFileName.value() != ::boost::none) {
    const std::string tag = params_.multiplyADFileName.value().value() +
      "_" + std::to_string(count_multiplyad_++);
    diagnostics(tag, fset);
  }

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void WriteVariances::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  if (params_.leftInverseFileName.value() != ::boost::none) {
    const std::string tag = params_.leftInverseFileName.value().value() +
      "_" + std::to_string(count_leftinversemultiply_++);
    diagnostics(tag, fset);
  }

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void WriteVariances::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
