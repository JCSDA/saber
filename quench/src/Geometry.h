/*
 * (C) Copyright 2022 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"

#include "eckit/mpi/Comm.h"

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

  // Vertical coordinate
  oops::OptionalParameter<std::vector<double>> vert_coord{"vert_coord", this};

  // Vertical coordinate from file
  oops::OptionalParameter<eckit::LocalConfiguration> vert_coordFromFile{"vert_coord from file",
    this};

  // Mask type
  oops::Parameter<std::string> maskType{"mask type", "none", this};

  // Mask path
  oops::Parameter<std::string> maskPath{"mask path", "../quench/data/landsea.nc", this};
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

  // Model data
  oops::Parameter<eckit::LocalConfiguration> modelData{"model data", eckit::LocalConfiguration(),
    this};

  // Aliases for model files
  oops::Parameter<std::vector<AliasParameters>> alias{"alias", {}, this};

  // Latitudes from south to north in files
  oops::Parameter<bool> latSouthToNorth{"latitude south to north", true, this};

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
  std::vector<size_t> variableSizes(const oops::Variables & vars) const;

  // Levels direction
  bool levelsAreTopDown() const
    {return levelsAreTopDown_;}

  // Accessors
  const eckit::mpi::Comm & getComm() const
    {return comm_;}
  const size_t halo() const
    {return halo_;}
  const atlas::Grid grid() const
    {return grid_;}
  const std::string gridType() const
    {return gridType_;}
  const atlas::grid::Partitioner partitioner() const
    {return partitioner_;}
  const atlas::Mesh mesh() const
    {return mesh_;}
  const atlas::FunctionSpace & functionSpace() const
    {return functionSpace_;}
  const atlas::FieldSet & fields() const
    {return fields_;}
  const size_t & levels(const size_t & groupIndex) const
    {return groups_[groupIndex].levels_;}
  const size_t & levels(const std::string & var) const
    {return groups_[groupIndex_.at(var)].levels_;}
  size_t groups() const
    {return groups_.size();}
  size_t groupIndex(const std::string & var) const
    {return groupIndex_.at(var);}
  const eckit::LocalConfiguration & modelData() const
    {return modelData_;}
  const std::vector<eckit::LocalConfiguration> & alias() const
    {return alias_;}
  const bool & latSouthToNorth() const
    {return latSouthToNorth_;}
  const eckit::LocalConfiguration & interpolation() const
    {return interpolation_;}
  const bool & duplicatePoints() const
    {return duplicatePoints_;}

 private:
  // Print
  void print(std::ostream &) const;

  // Read land-sea mask
  void readSeaMask(const std::string &,
                   const size_t &,
                   const std::string &,
                   atlas::Field &) const;

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
    size_t levels_;
    std::string lev2d_;
    atlas::Field vert_coord_;
    std::vector<double> vert_coord_avg_;
    double gmaskSize_;
  };

  // Geometry fields
  atlas::FieldSet fields_;

  // Groups
  std::vector<groupData> groups_;

  // Levels direction
  bool levelsAreTopDown_;

  // Model data configuration
  eckit::LocalConfiguration modelData_;

  // Aliases
  std::vector<eckit::LocalConfiguration> alias_;

  // Latitudes from south to north in files
  bool latSouthToNorth_;

  // Interpolation configuration
  eckit::LocalConfiguration interpolation_;

  // Duplicate points
  bool duplicatePoints_;
};

// -----------------------------------------------------------------------------

}  // namespace quench
