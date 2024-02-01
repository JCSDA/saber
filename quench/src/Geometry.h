/*
 * (C) Copyright 2022 UCAR.
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

#include "eckit/config/Configuration.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"

#include "eckit/mpi/Comm.h"

#include "oops/mpi/mpi.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

namespace quench {

// -----------------------------------------------------------------------------
/// Orography parameters

class OrographyParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OrographyParameters, Parameters)

 public:
  /// Top longitude [degrees]
  oops::RequiredParameter<double> topLon{"top longitude", this};

  /// Top latitude [degrees]
  oops::RequiredParameter<double> topLat{"top latitude", this};

  /// Zonal length [m]
  oops::RequiredParameter<double> zonalLength{"zonal length", this};

  /// Meridional length [m]
  oops::RequiredParameter<double> meridionalLength{"meridional length", this};

  /// Height (% of the bottom layer thickness, or absolute value if one level only)
  oops::RequiredParameter<double> height{"height", this};
};

// -----------------------------------------------------------------------------
/// Group parameters

class GroupParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GroupParameters, Parameters)

 public:
  /// Variables
  oops::RequiredParameter<std::vector<std::string>> variables{"variables", this};

  /// Number of levels
  oops::Parameter<size_t> levels{"levels", 1, this};

  /// Corresponding level for 2D variables (first or last)
  oops::Parameter<std::string> lev2d{"lev2d", "first", this};

  /// Orography
  oops::OptionalParameter<OrographyParameters> orography{"orography", this};

  /// Vertical coordinate
  oops::OptionalParameter<std::vector<double>> vert_coord{"vert_coord", this};

  /// Mask type
  oops::Parameter<std::string> maskType{"mask type", "none", this};

  /// Mask path
  oops::Parameter<std::string> maskPath{"mask path", "../quench/data/landsea.nc", this};
};

// -----------------------------------------------------------------------------
/// Geometry parameters

class GeometryParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeometryParameters, Parameters)

 public:
  /// Function space
  oops::RequiredParameter<std::string> functionSpace{"function space", this};

  /// Grid
  oops::RequiredParameter<eckit::LocalConfiguration> grid{"grid", this};

  /// Partitioner
  oops::Parameter<std::string> partitioner{"partitioner", "equal_regions", this};

  /// Variables groups
  oops::RequiredParameter<std::vector<GroupParameters>> groups{"groups", this};

  /// Halo size
  oops::Parameter<size_t> halo{"halo", 0, this};

  /// No point on last task
  oops::Parameter<bool> noPointOnLastTask{"no point on last task", false, this};
};

// -----------------------------------------------------------------------------
/// Geometry handles geometry for quench model.

class Geometry : public util::Printable,
                 private util::ObjectCounter<Geometry> {
 public:
  static const std::string classname() {return "quench::Geometry";}

  Geometry(const eckit::Configuration &,
           const eckit::mpi::Comm & comm = oops::mpi::world());
  Geometry(const Geometry &);

  const eckit::mpi::Comm & getComm() const {return comm_;}
  const size_t halo() const {return halo_;}
  const atlas::Grid grid() const {return grid_;}
  const std::string gridType() const {return gridType_;}
  const atlas::grid::Partitioner partitioner() const {return partitioner_;}
  const atlas::Mesh mesh() const {return mesh_;}
  const atlas::FunctionSpace & functionSpace() const {return functionSpace_;}
  const atlas::FieldSet & fields() const {return groups_[0].fields_;}
  const atlas::FieldSet & fields(const size_t & groupIndex) const
    {return groups_[groupIndex].fields_;}
  size_t levels(const size_t & groupIndex) const {return groups_[groupIndex].levels_;}
  size_t levels(const std::string & var) const;
  size_t groups() const {return groups_.size();}
  size_t groupIndex(const std::string & var) const;

  size_t variableSize(const std::string &) const;
  size_t maskLevel(const std::string &, const size_t &) const;
  std::vector<size_t> variableSizes(const oops::Variables & vars) const;
  void latlon(std::vector<double> &, std::vector<double> &, const bool) const;
  bool levelsAreTopDown() const {return true;}

 private:
  void print(std::ostream &) const;
  void readSeaMask(const std::string &, const size_t &, const std::string &, atlas::Field &) const;
  const eckit::mpi::Comm & comm_;
  size_t halo_;
  atlas::Grid grid_;
  std::string gridType_;
  bool regionalGrid_;
  atlas::grid::Partitioner partitioner_;
  atlas::Mesh mesh_;
  atlas::FunctionSpace functionSpace_;
  std::unordered_map<std::string, size_t> groupIndex_;
  struct groupData {
    size_t levels_;
    std::string lev2d_;
    std::vector<double> vert_coord_;
    atlas::FieldSet fields_;
    double gmaskSize_;
  };
  std::vector<groupData> groups_;
};
// -----------------------------------------------------------------------------

}  // namespace quench
