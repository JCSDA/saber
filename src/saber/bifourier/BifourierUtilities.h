/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variable.h"

#include "saber/bifourier/BifourierTransformBase.h"

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

void createField1D(const std::string & prefix,
                   const oops::Variable &,
                   atlas::FieldSet &);

// -----------------------------------------------------------------------------

void createField2D(const std::string & prefix,
                   const size_t &,
                   const oops::Variable &,
                   atlas::FieldSet &);

// -----------------------------------------------------------------------------

void createField2D(const std::string & prefix,
                   const oops::Variable &,
                   atlas::FieldSet &);

// -----------------------------------------------------------------------------

void createField2D(const std::string & prefix,
                   const oops::Variable &,
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

atlas::array::ArrayView<double, 1> getView1D(const std::string &,
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

atlas::array::ArrayView<double, 2> getView2D(const std::string &,
                                             const oops::Variable &,
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

void defNcVar(const std::string &,
              const oops::Variable &,
              const int &,
              const std::vector<int> &,
              std::vector<int> &,
              size_t &);

// -----------------------------------------------------------------------------

void defNcVar(const std::string &,
              const oops::Variable &,
              const oops::Variable &,
              const int &,
              const std::vector<int> &,
              std::vector<int> &,
              size_t &);

// -----------------------------------------------------------------------------

void putNcVar(const std::string &,
              const oops::Variable &,
              const int &,
              const std::vector<int> &,
              size_t &,
              const atlas::FieldSet &,
              const BifourierTransformBase &);

// -----------------------------------------------------------------------------

void putNcVar(const std::string &,
              const oops::Variable &,
              const oops::Variable &,
              const int &,
              const std::vector<int> &,
              size_t &,
              const atlas::FieldSet &,
              const BifourierTransformBase &);

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
