/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_UTIL_ATLASFIELDSET_H_
#define SABER_UTIL_ATLASFIELDSET_H_

#include "atlas/field.h"

#include "oops/util/Logger.h"

namespace saber {

// -----------------------------------------------------------------------------

void multiplyAtlasFieldSet(atlas::FieldSet & fset, const atlas::FieldSet & mulFset) {
  oops::Log::trace() << "multiplyAtlasFieldSet starting" << std::endl;
  // Loop over multiplier fields
  for (auto mulField : mulFset) {
    // Get increment field with the same name
    atlas::Field field = fset.field(mulField.name());

    // Get data and multiply
    if (mulField.rank() == 1) {
      const auto mulView = atlas::array::make_view<double, 1>(mulField);
      auto view = atlas::array::make_view<double, 1>(field);
      for (int j0 = 0; j0 < mulField.shape(0); ++j0) {
        view(j0) = view(j0)*mulView(j0);
      }
    }
    if (mulField.rank() == 2) {
      const auto mulView = atlas::array::make_view<double, 2>(mulField);
      auto view = atlas::array::make_view<double, 2>(field);
      for (int j0 = 0; j0 < mulField.shape(0); ++j0) {
        for (int j1 = 0; j1 < mulField.shape(1); ++j1) {
          view(j0, j1) = view(j0, j1)*mulView(j0, j1);
        }
      }
    }
    if (mulField.rank() == 3) {
      const auto mulView = atlas::array::make_view<double, 3>(mulField);
      auto view = atlas::array::make_view<double, 3>(field);
      for (int j0 = 0; j0 < mulField.shape(0); ++j0) {
        for (int j1 = 0; j1 < mulField.shape(1); ++j1) {
          for (int j2 = 0; j2 < mulField.shape(2); ++j2) {
            view(j0, j1, j2) = view(j0, j1, j2)*mulView(j0, j1, j2);
          }
        }
      }
    }
  }
  oops::Log::trace() << "multiplyAtlasFieldSet starting" << std::endl;
}

// -----------------------------------------------------------------------------

void divideAtlasFieldSet(atlas::FieldSet & fset, const  atlas::FieldSet & divFset) {
  oops::Log::trace() << "divideAtlasFieldSet starting" << std::endl;
  // Loop over divider fields
  for (auto divField : divFset) {
    // Get increment field with the same name
    atlas::Field field = fset.field(divField.name());

    // Get data and divide
    if (divField.rank() == 1) {
      const auto divView = atlas::array::make_view<double, 1>(divField);
      auto view = atlas::array::make_view<double, 1>(field);
      for (int j0 = 0; j0 < divField.shape(0); ++j0) {
        view(j0) = view(j0)/divView(j0);
      }
    }
    if (divField.rank() == 2) {
      const auto divView = atlas::array::make_view<double, 2>(divField);
      auto view = atlas::array::make_view<double, 2>(field);
      for (int j0 = 0; j0 < divField.shape(0); ++j0) {
        for (int j1 = 0; j1 < divField.shape(1); ++j1) {
          view(j0, j1) = view(j0, j1)/divView(j0, j1);
        }
      }
    }
    if (divField.rank() == 3) {
      const auto divView = atlas::array::make_view<double, 3>(divField);
      auto view = atlas::array::make_view<double, 3>(field);
      for (int j0 = 0; j0 < divField.shape(0); ++j0) {
        for (int j1 = 0; j1 < divField.shape(1); ++j1) {
          for (int j2 = 0; j2 < divField.shape(2); ++j2) {
            view(j0, j1, j2) = view(j0, j1, j2)/divView(j0, j1, j2);
          }
        }
      }
    }
  }
  oops::Log::trace() << "divideAtlasFieldSet starting" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_UTIL_ATLASFIELDSET_H_
