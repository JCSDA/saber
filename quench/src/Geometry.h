/*
 * (C) Copyright 2022 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"

#include "eckit/mpi/Comm.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace quench {

// -----------------------------------------------------------------------------
/// Orography parameters

class OrographyParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OrographyParameters, Parameters)

 public:
  // Top longitude [degrees]
  oops::RequiredParameter<double> topLon{"top longitude", this};

  // Top latitude [degrees]
  oops::RequiredParameter<double> topLat{"top latitude", this};

  // Zonal length [m]
  oops::RequiredParameter<double> zonalLength{"zonal length", this};

  // Meridional length [m]
  oops::RequiredParameter<double> meridionalLength{"meridional length", this};

  // Height (% of the bottom layer thickness, or absolute value if one level only)
  oops::RequiredParameter<double> height{"height", this};
};

// -----------------------------------------------------------------------------
/// Group parameters

class GroupParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GroupParameters, Parameters)

 public:
  // Variables
  oops::RequiredParameter<std::vector<std::string>> variables{"variables", this};

  // Number of levels
  oops::Parameter<size_t> levels{"levels", 1, this};

  // Corresponding level for 2D variables (first or last)
  oops::Parameter<std::string> lev2d{"lev2d", "first", this};

  // Orography
  oops::OptionalParameter<OrographyParameters> orography{"orography", this};

  // Vertical coordinate configuration
  oops::OptionalParameter<eckit::LocalConfiguration> vertCoordConf{
    "vertical coordinate", this};

  // Mask type
  oops::Parameter<std::string> maskType{"mask type", "none", this};

  // Mask path
  oops::OptionalParameter<std::string> maskPath{"mask path", this};
};

// -----------------------------------------------------------------------------
/// Alias elemental paramaters

class AliasParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(AliasParameters, oops::Parameters)

 public:
  // In code
  oops::RequiredParameter<std::string> inCode{"in code", this};
  // In model file
  oops::RequiredParameter<std::string> inFile{"in file", this};
  // Optional parameters for States transformations
  // Scaling factor (e.g. for units conversion)
  oops::OptionalParameter<double> scalingFactor{"scaling factor", this};
  // Toggle log10 transformation
  oops::OptionalParameter<bool> logTransf{"log transform", this};
  // Additive constant (prior to log10 transformation)
  oops::OptionalParameter<double> addConst{"additive constant", this};
};

// -----------------------------------------------------------------------------
/// Interpolation paramaters

class InterpolationParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(InterpolationParameters, oops::Parameters)

 public:
  // Interpolation type
  oops::RequiredParameter<std::string> interpType{"interpolation type", this};
};

// -----------------------------------------------------------------------------
/// Geometry parameters

class GeometryParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeometryParameters, Parameters)

 public:
  // Function space
  oops::RequiredParameter<std::string> functionSpace{"function space", this};

  // Grid
  oops::RequiredParameter<eckit::LocalConfiguration> grid{"grid", this};

  // Partitioner
  oops::Parameter<std::string> partitioner{"partitioner", "equal_regions", this};

  // Variables groups
  oops::RequiredParameter<std::vector<GroupParameters>> groups{"groups", this};

  // Halo size
  oops::Parameter<size_t> halo{"halo", 0, this};

  // No point on last task
  oops::Parameter<bool> noPointOnLastTask{"no point on last task", false, this};

  // Levels top-down
  oops::Parameter<bool> levelsAreTopDown{"levels are top down", true, this};

  // Levels counter origin
  oops::Parameter<size_t> levelsCountFrom{"levels count from", 1, this};

  // Model data
  oops::Parameter<eckit::LocalConfiguration> modelData{"model data", eckit::LocalConfiguration(),
    this};

  // Variables name alias for model files
  oops::Parameter<std::vector<AliasParameters>> alias{"alias", {}, this};

  // Check longitudes/latitudes from file
  oops::OptionalParameter<eckit::LocalConfiguration> checkLonLat{"check lon/lat from file", this};

  // IO parameters
  oops::Parameter<eckit::LocalConfiguration> io{"io", eckit::LocalConfiguration(), this};

  // Interpolation parameters
  oops::OptionalParameter<InterpolationParameters> interpolation{"interpolation", this};
};

// -----------------------------------------------------------------------------
/// Geometry class

class Geometry : public util::Printable,
                 private util::ObjectCounter<Geometry> {
 public:
  static const std::string classname()
    {return "quench::Geometry";}

  // Constructors
  Geometry(const eckit::Configuration &,
           const eckit::mpi::Comm & comm = oops::mpi::world());
  Geometry(const Geometry &);

  // Variables sizes
  std::vector<size_t> variableSizes(const oops::Variables &) const;
  std::vector<size_t> variableSizes(const std::vector<std::string> &) const;

  // Levels direction
  const bool levelsAreTopDown() const
    {return levelsAreTopDown_;}

  // Levels counter origin
  const bool levelsCountFrom() const
    {return levelsCountFrom_;}

  // Accessors
  const eckit::mpi::Comm & getComm() const
    {return comm_;}
  size_t halo() const
    {return halo_;}
  atlas::Grid grid() const
    {return grid_;}
  std::string gridType() const
    {return gridType_;}
  atlas::grid::Partitioner partitioner() const
    {return partitioner_;}
  atlas::Mesh mesh() const
    {return mesh_;}
  const atlas::FunctionSpace & functionSpace() const
    {return functionSpace_;}
  const atlas::FieldSet & fields() const
    {return fields_;}
  size_t levels(const size_t & groupIndex) const
    {return groups_[groupIndex].levels_;}
  size_t levels(const std::string & var) const
    {return groups_[groupIndex(var)].levels_;}
  size_t groups() const
    {return groups_.size();}
  size_t groupIndex(const std::string & var) const;
  const eckit::LocalConfiguration & modelData() const
    {return modelData_;}
  const std::vector<eckit::LocalConfiguration> & alias() const
    {return alias_;}
  const eckit::LocalConfiguration & io() const
    {return io_;}
  const eckit::LocalConfiguration & interpolation() const
    {return interpolation_;}
  bool duplicatePoints() const
    {return duplicatePoints_;}
  const std::vector<double> & vertCoordAvg(const std::string & var) const
    {return groups_[groupIndex(var)].vertCoordAvg_;}
  const oops::GeometryData & generic() const
    {return *geomData_;}

 private:
  // Communicator
  const eckit::mpi::Comm & comm_;

  // Halo size
  size_t halo_;

  // ATLAS grid
  atlas::Grid grid_;

  // ATLAS grid type
  std::string gridType_;

  // ATLAS grid partitioner
  atlas::grid::Partitioner partitioner_;

  // ATLAS mesh
  atlas::Mesh mesh_;

  // ATLAS function space
  atlas::FunctionSpace functionSpace_;

  // Group name to group index mapping
  std::unordered_map<std::string, size_t> groupIndex_;

  // Group data structure
  struct groupData {
    GroupParameters params_;
    size_t index_;
    size_t levels_;
    std::string lev2d_;
    atlas::Field vertCoord_;
    std::vector<double> vertCoordAvg_;
    double gmaskSize_;
  };

  // Geometry fields
  atlas::FieldSet fields_;

  // Groups
  std::vector<groupData> groups_;

  // Levels direction
  bool levelsAreTopDown_;

  // Levels counter origin
  size_t levelsCountFrom_;

  // Model data configuration
  eckit::LocalConfiguration modelData_;

  // Variables name alias
  std::vector<eckit::LocalConfiguration> alias_;

  // IO configuration
  eckit::LocalConfiguration io_;

  // Interpolation configuration
  eckit::LocalConfiguration interpolation_;

  // Duplicate points
  bool duplicatePoints_;

  // Geometry data structure
  std::unique_ptr<oops::GeometryData> geomData_;

  // Private methods

  // Print
  void print(std::ostream &) const;

  // Setup alias
  void setupAlias(const GeometryParameters &);

  // Setup group vertical coordinate
  void setupVertCoord(groupData &);

  // Setup group mask
  void setupMask(groupData &);

  // Check longitudes/latitudes from file
  void checkLonLat(const eckit::Configuration &);
};

// -----------------------------------------------------------------------------

}  // namespace quench
