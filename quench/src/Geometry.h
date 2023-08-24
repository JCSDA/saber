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

  /// Vertical unit
  oops::OptionalParameter<std::vector<double>> vunit{"vunit", this};

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
  oops::OptionalParameter<size_t> halo{"halo", this};

  /// No point on last task
  oops::Parameter<bool> noPointOnLastTask{"no point on last task", false, this};
};

// -----------------------------------------------------------------------------
/// Geometry handles geometry for  model.

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
  atlas::FunctionSpace & functionSpace() {return functionSpace_;}
  const atlas::FieldSet & extraFields() const {return groups_[0].extraFields_;}
  atlas::FieldSet & extraFields() {return groups_[0].extraFields_;}
  const atlas::FieldSet & extraFields(const size_t & groupIndex) const
    {return groups_[groupIndex].extraFields_;}
  size_t levels(const size_t & groupIndex) const {return groups_[groupIndex].levels_;}
  size_t levels(const std::string & var) const;
  size_t groups() const {return groups_.size();}
  size_t groupIndex(const std::string & var) const;

  size_t variableSize(const std::string &) const;
  size_t maskLevel(const std::string &, const size_t &) const;
  std::vector<size_t> variableSizes(const oops::Variables & vars) const;
  void latlon(std::vector<double> &, std::vector<double> &, const bool) const;
  bool levelsAreTopDown() const {return true;}
  bool iodaBased() const {return iodaBased_;}

 private:
  void print(std::ostream &) const;
  void readSeaMask(const std::string &, const size_t &, const std::string &, atlas::Field &) const;
  const eckit::mpi::Comm & comm_;
  size_t halo_;
  atlas::Grid grid_;
  std::string gridType_;
  bool unstructuredGrid_;
  atlas::grid::Partitioner partitioner_;
  atlas::grid::Distribution distribution_;
  atlas::Mesh mesh_;
  atlas::FunctionSpace functionSpace_;
  std::unordered_map<std::string, size_t> groupIndex_;
  struct groupData {
    size_t levels_;
    std::string lev2d_;
    std::vector<double> vunit_;
    atlas::FieldSet extraFields_;
    double gmaskSize_;
  };
  std::vector<groupData> groups_;
  bool iodaBased_;
};
// -----------------------------------------------------------------------------

}  // namespace quench
