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

void multiplyAtlasFieldSet(atlas::FieldSet * fieldSet, atlas::FieldSet * multiplier) {
  oops::Log::trace() << "multiplyAtlasFieldSet starting" << std::endl;
  // Loop over multiplier fields
  for (atlas::FieldSet::iterator it = multiplier->begin(); it != multiplier->end(); ++it) {
    // Get increment field with the same name
    atlas::Field field = fieldSet->field(it->name());

    // Get data and multiply
    if (it->rank() == 1) {
      const auto multiplierView = atlas::array::make_view<double, 1>(*it);
      auto view = atlas::array::make_view<double, 1>(field);
      for (int j0 = 0; j0 < it->shape(0); ++j0) {
        view(j0) = view(j0)*multiplierView(j0);
      }
    }
    if (it->rank() == 2) {
      const auto multiplierView = atlas::array::make_view<double, 2>(*it);
      auto view = atlas::array::make_view<double, 2>(field);
      for (int j0 = 0; j0 < it->shape(0); ++j0) {
        for (int j1 = 0; j1 < it->shape(1); ++j1) {
          view(j0, j1) = view(j0, j1)*multiplierView(j0, j1);
        }
      }
    }
    if (it->rank() == 3) {
      const auto multiplierView = atlas::array::make_view<double, 3>(*it);
      auto view = atlas::array::make_view<double, 3>(field);
      for (int j0 = 0; j0 < it->shape(0); ++j0) {
        for (int j1 = 0; j1 < it->shape(1); ++j1) {
          for (int j2 = 0; j2 < it->shape(2); ++j2) {
            view(j0, j1, j2) = view(j0, j1, j2)*multiplierView(j0, j1, j2);
          }
        }
      }
    }
  }
  oops::Log::trace() << "multiplyAtlasFieldSet starting" << std::endl;
}

// -----------------------------------------------------------------------------

void divideAtlasFieldSet(atlas::FieldSet * fieldSet, atlas::FieldSet * divider) {
  oops::Log::trace() << "divideAtlasFieldSet starting" << std::endl;
  // Loop over divider fields
  for (atlas::FieldSet::iterator it = divider->begin(); it != divider->end(); ++it) {
    // Get increment field with the same name
    atlas::Field field = fieldSet->field(it->name());

    // Get data and divide
    if (it->rank() == 1) {
      const auto dividerView = atlas::array::make_view<double, 1>(*it);
      auto view = atlas::array::make_view<double, 1>(field);
      for (int j0 = 0; j0 < it->shape(0); ++j0) {
        view(j0) = view(j0)/dividerView(j0);
      }
    }
    if (it->rank() == 2) {
      const auto dividerView = atlas::array::make_view<double, 2>(*it);
      auto view = atlas::array::make_view<double, 2>(field);
      for (int j0 = 0; j0 < it->shape(0); ++j0) {
        for (int j1 = 0; j1 < it->shape(1); ++j1) {
          view(j0, j1) = view(j0, j1)/dividerView(j0, j1);
        }
      }
    }
    if (it->rank() == 3) {
      const auto dividerView = atlas::array::make_view<double, 3>(*it);
      auto view = atlas::array::make_view<double, 3>(field);
      for (int j0 = 0; j0 < it->shape(0); ++j0) {
        for (int j1 = 0; j1 < it->shape(1); ++j1) {
          for (int j2 = 0; j2 < it->shape(2); ++j2) {
            view(j0, j1, j2) = view(j0, j1, j2)/dividerView(j0, j1, j2);
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
