/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Variables.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace saber {

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
std::vector<atlas::FieldSet> readInputFields(
  const oops::Geometry<MODEL> & resol,
  const oops::Variables & vars,
  const util::DateTime & date,
  const boost::optional<std::vector<eckit::LocalConfiguration>> & inputFields) {
  // Vector of FieldSets
  std::vector<atlas::FieldSet> fsetVec;

  // Get number of MPI tasks and OpenMP threads
  std::string mpi(std::to_string(resol.getComm().size()));
  std::string omp("1");
  # pragma omp parallel
  {
    omp = std::to_string(omp_get_num_threads());
  }

  // Read block input fields function
  if (inputFields != boost::none) {
    if (inputFields->size() > 0) {
      // Create Increment
      oops::Increment<MODEL> dx(resol, vars, date);

      // Loop over block input fields
      for (const auto & inputField : *inputFields) {
        // Get input field file configuration
        eckit::LocalConfiguration file = inputField.getSubConfiguration("file");

        // Replace patterns
        util::seekAndReplace(file, "_MPI_", mpi);
        util::seekAndReplace(file, "_OMP_", omp);

        // Read block input field as Increment
        dx.read(file);

        // Define FieldSet name
        std::string name = inputField.getString("parameter");
        if (inputField.has("component")) {
          name += "::" + std::to_string(inputField.getInt("component"));
        }

        // Transform Increment into FieldSet
        oops::Log::test() << "Norm of input parameter " << name << ": " << dx.norm() << std::endl;
        atlas::FieldSet fset;
        fset.name() = name;
        for (const atlas::Field field : dx.fieldSet()) {
          fset.add(field);
        }
        fsetVec.push_back(fset);
      }
    }
  }
  return fsetVec;
}

// -----------------------------------------------------------------------------

}  // namespace saber
