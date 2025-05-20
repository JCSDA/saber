/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierUtilities.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/utils/MD5.h"

#include "oops/util/Logger.h"

using atlas::array::make_datatype;
using atlas::array::make_shape;
using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

std::string generateSpectralUid(const atlas::FunctionSpace & fspace,
                                const eckit::mpi::Comm & comm) {
  oops::Log::trace() << "saber::bifourier::generateSpectralUid starting" << std::endl;

  // Check function space type
  ASSERT(fspace.type() == "PointCloud");

  // Create separate lon/lat field (because of a bug in ATLAS)
  atlas::FieldSet fset;
  atlas::Field lonField = fspace.createField<double>(atlas::option::name("lon"));
  atlas::Field latField = fspace.createField<double>(atlas::option::name("lat"));
  fset.add(lonField);
  fset.add(latField);
  auto lonView = make_view<double, 1>(lonField);
  auto latView = make_view<double, 1>(latField);
  const auto lonlatView = make_view<double, 2>(fspace.lonlat());
  for (int jj = 0; jj < lonField.shape(0); ++jj) {
    lonView(jj) = lonlatView(jj, 0);
    latView(jj) = lonlatView(jj, 1);
  }

  // Prepare global lon/lat fields
  atlas::FieldSet fsetGlb;
  atlas::Field lonGlbField = fspace.createField<double>(atlas::option::name("lon") |
    atlas::option::global());
  atlas::Field latGlbField = fspace.createField<double>(atlas::option::name("lat") |
    atlas::option::global());
  fsetGlb.add(lonGlbField);
  fsetGlb.add(latGlbField);

  // Gather lon/lat
  fspace.gather(fset, fsetGlb);

  // Compute UID
  std::string uid;
  size_t uidLength;
  if (comm.rank() == 0) {
    const auto lonGlbView = make_view<double, 1>(lonGlbField);
    const auto latGlbView = make_view<double, 1>(latGlbField);
    eckit::MD5 hash;
    const double hashTrunc = 1.0e8;
    for (int jj = 0; jj < lonGlbField.shape(0); ++jj) {
      hash.add(std::round(lonGlbView(jj)*hashTrunc));
      hash.add(std::round(latGlbView(jj)*hashTrunc));
    }
    uid = hash.digest();
    uidLength = uid.size();
  }

  // Broadcast UID length
  comm.broadcast(uidLength, 0);

  // Broadcast UID
  if (comm.rank() != 0) {
    uid.resize(uidLength);
  }
  comm.broadcast(uid.begin(), uid.end(), 0);

  oops::Log::trace() << "saber::bifourier::generateSpectralUid done" << std::endl;
  return uid;
}

// -----------------------------------------------------------------------------

std::string fieldName(const std::string & prefix,
                      const oops::Variable & var) {
  oops::Log::trace() << "saber::bifourier::fieldName starting" << std::endl;

  const std::string name = prefix + "_" + var.name();

  oops::Log::trace() << "saber::bifourier::fieldName done" << std::endl;
  return name;
}

// -----------------------------------------------------------------------------

std::string fieldName(const std::string & prefix,
                      const oops::Variable & outputVar,
                      const oops::Variable & inputVar) {
  oops::Log::trace() << "saber::bifourier::fieldName starting" << std::endl;

  const std::string name = prefix + "_" + outputVar.name() + "_from_" + inputVar.name();

  oops::Log::trace() << "saber::bifourier::fieldName done" << std::endl;
  return name;
}

// -----------------------------------------------------------------------------

void createFieldProfile(const std::string & prefix,
                        const oops::Variable & var,
                        atlas::FieldSet & fset) {
  oops::Log::trace() << "saber::bifourier::createFieldProfile starting" << std::endl;

  // Get matrix field name
  const std::string name = fieldName(prefix, var);

  // Create matrix field
  atlas::Field field(name, make_datatype<double>(),
    make_shape(var.getLevels()));

  // Get matrix view
  auto view = make_view<double, 1>(field);

  // Set matrix to zero
  view.assign(0.0);

  // Add matrix field
  fset.add(field);

  oops::Log::trace() << "saber::bifourier::createFieldProfile done" << std::endl;
}

// -----------------------------------------------------------------------------

void createField2D(const std::string & prefix,
                   const size_t & bins,
                   const oops::Variable & var,
                   atlas::FieldSet & fset) {
  oops::Log::trace() << "saber::bifourier::createField2D starting" << std::endl;

  // Get field name
  const std::string name = fieldName(prefix, var);

  // Create field
  atlas::Field field(name, make_datatype<double>(), make_shape(bins, var.getLevels()));

  // Get output perturbation view
  auto view = make_view<double, 2>(field);

  // Set field to zero
  view.assign(0.0);

  // Add field
  fset.add(field);

  oops::Log::trace() << "saber::bifourier::createField2D done" << std::endl;
}

// -----------------------------------------------------------------------------

void createField3D(const std::string & prefix,
                   const size_t & bins,
                   const oops::Variable & var,
                   atlas::FieldSet & fset) {
  oops::Log::trace() << "saber::bifourier::createField3D starting" << std::endl;

  // Get matrix field name
  const std::string name = fieldName(prefix, var);

  // Create matrix field
  atlas::Field field(name, make_datatype<double>(),
    make_shape(bins, var.getLevels(), var.getLevels()));

  // Get matrix view
  auto view = make_view<double, 3>(field);

  // Set matrix to zero
  view.assign(0.0);

  // Add matrix field
  fset.add(field);

  oops::Log::trace() << "saber::bifourier::createField3D done" << std::endl;
}

// -----------------------------------------------------------------------------

void createField3D(const std::string & prefix,
                   const size_t & bins,
                   const oops::Variable & outputVar,
                   const oops::Variable & inputVar,
                   atlas::FieldSet & fset) {
  oops::Log::trace() << "saber::bifourier::createField3D starting" << std::endl;

  // Get matrix field name
  const std::string name = fieldName(prefix, outputVar, inputVar);

  // Create matrix field
  atlas::Field field(name, make_datatype<double>(),
    make_shape(bins, outputVar.getLevels(), inputVar.getLevels()));

  // Get matrix view
  auto view = make_view<double, 3>(field);

  // Set matrix to zero
  view.assign(0.0);

  // Add matrix field
  fset.add(field);

  oops::Log::trace() << "saber::bifourier::createField3D done" << std::endl;
}

// -----------------------------------------------------------------------------

atlas::Field getField(const std::string & prefix,
                      const oops::Variable & var,
                      const atlas::FieldSet & fset) {
  oops::Log::trace() << "saber::bifourier::getField starting" << std::endl;

  // Get matrix field name
  const std::string name = fieldName(prefix, var);

  // Get matrix field
  auto field = fset[name];

  // Check number of levels
  if (field.rank() == 1) {
    // Profile
    ASSERT(field.shape(0) == static_cast<int>(var.getLevels()));
  } else if (field.rank() == 2) {
    // 2D field
    ASSERT(field.shape(1) == static_cast<int>(var.getLevels()));
  } else if (field.rank() == 3) {
    // 3D field
    ASSERT(field.shape(1) == static_cast<int>(var.getLevels()));
    ASSERT(field.shape(2) == static_cast<int>(var.getLevels()));
  } else {
    // Wrong rank
    throw eckit::Exception("wrong rank", Here());
  }

  oops::Log::trace() << "saber::bifourier::getField done" << std::endl;
  return field;
}

// -----------------------------------------------------------------------------

atlas::Field getField(const std::string & prefix,
                      const oops::Variable & outputVar,
                      const oops::Variable & inputVar,
                      const atlas::FieldSet & fset) {
  oops::Log::trace() << "saber::bifourier::getField starting" << std::endl;

  // Get matrix field name
  const std::string name = fieldName(prefix, outputVar, inputVar);

  // Get matrix field
  auto field = fset[name];

  // Check number of levels
  ASSERT(field.shape(1) == static_cast<int>(outputVar.getLevels()));
  ASSERT(field.shape(2) == static_cast<int>(inputVar.getLevels()));

  oops::Log::trace() << "saber::bifourier::getField done" << std::endl;
  return field;
}

// -----------------------------------------------------------------------------

atlas::array::ArrayView<double, 1> getViewProfile(const std::string & prefix,
                                                  const oops::Variable & var,
                                                  const atlas::FieldSet & fset) {
  oops::Log::trace() << "saber::bifourier::getViewProfile starting" << std::endl;

  // Get matrix field
  auto field = getField(prefix, var, fset);

  // Get matrix view
  auto view = make_view<double, 1>(field);

  oops::Log::trace() << "saber::bifourier::getViewProfile done" << std::endl;
  return view;
}

// -----------------------------------------------------------------------------

atlas::array::ArrayView<double, 2> getView2D(const oops::Variable & var,
                                             const oops::FieldSet3D & fset) {
  oops::Log::trace() << "saber::bifourier::getView2D starting" << std::endl;

  // Get field name
  auto field = fset[var.name()];;

  // Get view
  auto view = make_view<double, 2>(field);

  oops::Log::trace() << "saber::bifourier::getView2D done" << std::endl;
  return view;
}

// -----------------------------------------------------------------------------

atlas::array::ArrayView<double, 2> getView2D(const std::string & prefix,
                                             const oops::Variable & var,
                                             const atlas::FieldSet & fset) {
  oops::Log::trace() << "saber::bifourier::getView2D starting" << std::endl;

  // Get field name
  auto field = getField(prefix, var, fset);

  // Get view
  auto view = make_view<double, 2>(field);

  oops::Log::trace() << "saber::bifourier::getView2D done" << std::endl;
  return view;
}

// -----------------------------------------------------------------------------

atlas::array::ArrayView<double, 3> getView3D(const std::string & prefix,
                                             const oops::Variable & var,
                                             const atlas::FieldSet & fset) {
  oops::Log::trace() << "saber::bifourier::getView3D starting" << std::endl;

  // Get matrix field
  auto field = getField(prefix, var, fset);

  // Get matrix view
  auto view = make_view<double, 3>(field);

  oops::Log::trace() << "saber::bifourier::getView3D done" << std::endl;
  return view;
}

// -----------------------------------------------------------------------------

atlas::array::ArrayView<double, 3> getView3D(const std::string & prefix,
                                             const oops::Variable & outputVar,
                                             const oops::Variable & inputVar,
                                             const atlas::FieldSet & fset) {
  oops::Log::trace() << "saber::bifourier::getView3D starting" << std::endl;

  // Get matrix field
  auto field = getField(prefix, outputVar, inputVar, fset);

  // Get matrix view
  auto view = make_view<double, 3>(field);

  oops::Log::trace() << "saber::bifourier::getView3D done" << std::endl;
  return view;
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
