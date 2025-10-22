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
#include "oops/mpi/mpi.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "src/GeometryParameters.h"
#include "src/Interpolation.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace quench {

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
  size_t groupIndex(const std::string &) const;
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

  // Interpolation
  Interpolation & getInterpolation(const Geometry &) const;

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
  mutable std::unique_ptr<oops::GeometryData> geomData_;

  // Interpolations vector
  mutable std::unordered_map<std::string, std::shared_ptr<Interpolation>> interpolations_;

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
