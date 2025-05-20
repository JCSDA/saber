/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variable.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

std::string generateSpectralUid(const atlas::FunctionSpace &,
                                const eckit::mpi::Comm &);

// -----------------------------------------------------------------------------

std::string fieldName(const std::string &,
                      const oops::Variable &);

// -----------------------------------------------------------------------------

std::string fieldName(const std::string &,
                      const oops::Variable &,
                      const oops::Variable &);

// -----------------------------------------------------------------------------

void createFieldProfile(const std::string & prefix,
                        const oops::Variable &,
                        atlas::FieldSet &);

// -----------------------------------------------------------------------------

void createField2D(const std::string & prefix,
                   const size_t &,
                   const oops::Variable &,
                   atlas::FieldSet &);

// -----------------------------------------------------------------------------

void createField3D(const std::string &,
                   const size_t &,
                   const oops::Variable &,
                   atlas::FieldSet &);

// -----------------------------------------------------------------------------

void createField3D(const std::string &,
                   const size_t &,
                   const oops::Variable &,
                   const oops::Variable &,
                   atlas::FieldSet &);

// -----------------------------------------------------------------------------

atlas::Field getField(const std::string &,
                      const oops::Variable &,
                      const atlas::FieldSet &);

// -----------------------------------------------------------------------------

atlas::Field getField(const std::string &,
                      const oops::Variable &,
                      const oops::Variable &,
                      const atlas::FieldSet &);

// -----------------------------------------------------------------------------

atlas::array::ArrayView<double, 1> getViewProfile(const std::string &,
                                                  const oops::Variable &,
                                                  const atlas::FieldSet &);

// -----------------------------------------------------------------------------

atlas::array::ArrayView<double, 2> getView2D(const oops::Variable &,
                                             const oops::FieldSet3D &);

// -----------------------------------------------------------------------------

atlas::array::ArrayView<double, 2> getView2D(const std::string &,
                                             const oops::Variable &,
                                             const atlas::FieldSet &);

// -----------------------------------------------------------------------------

atlas::array::ArrayView<double, 3> getView3D(const std::string &,
                                             const oops::Variable &,
                                             const atlas::FieldSet &);

// -----------------------------------------------------------------------------

atlas::array::ArrayView<double, 3> getView3D(const std::string &,
                                             const oops::Variable &,
                                             const oops::Variable &,
                                             const atlas::FieldSet &);

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
