/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/generic/OrographicInterp.h"

#include <algorithm>
#include <limits>

#include "eckit/exception/Exceptions.h"

using atlas::array::make_view;

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<OrographicInterp> makerOrographicInterp_("OrographicInterp");

// -----------------------------------------------------------------------------

OrographicInterp::OrographicInterp(const oops::GeometryData & outerGeometryData,
                                   const oops::Variables & outerVars,
                                   const eckit::Configuration & covarConfig,
                                   const Parameters_ & params,
                                   const oops::FieldSet3D & xb,
                                   const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    comm_(outerGeometryData.comm()),
    innerVars_(outerVars),
    params_(params)
{
  oops::Log::trace() << classname() << "::OrographicInterp starting" << std::endl;

  // Get ghost view
  const auto ghostView = make_view<int, 1>(outerGeometryData.functionSpace().ghost());

  // Vertical coordinates fieldset
  atlas::FieldSet vcFset;

  // List what vertical coordinates will be used
  for (const auto & var : innerVars_) {
    if (params_.fieldsMetaData.value().has(var.name())) {
      // Add interpolated variable
      interpVars_.push_back(var);

      // Get vertical coordinate name
      const std::string vcName = params_.fieldsMetaData.value().getString(
        var.name() + ".vert_coord");

      if (!vcFset.has(vcName)) {
        if (xb.has(vcName)) {
          // From background
          vcFset.add(xb[vcName]);
        } else {
          // From geometry fieldset
          if (outerGeometryData.fieldSet().has(vcName)) {
            vcFset.add(outerGeometryData.fieldSet()[vcName]);
          } else {
            throw eckit::UserError("variable " + vcName + " missing in geometry fields", Here());
          }
        }
      }
    }
  }

  // Process vertical coordinates
  for (const auto & vcField : vcFset) {
    oops::Log::info() << "Info     : Process vertical coordinate " << vcField.name() << std::endl;

    // Get number of levels
    const size_t nz = vcField.shape(1);
    ASSERT(nz > 1);

    // Get vertical coordinate field view
    const auto vcView = make_view<double, 2>(vcField);

    // Compute vertical coordinates min/max/avg profile
    std::vector<double> vcMin(nz, std::numeric_limits<double>().max());
    std::vector<double> vcMax(nz, -std::numeric_limits<double>().max());
    std::vector<double> vcAvg(nz, 0.0);
    int vcSize = 0.0;
    for (int jnode = 0; jnode < vcField.shape(0); ++jnode) {
      if (ghostView(jnode) == 0) {
        for (size_t jz = 0; jz < nz; ++jz) {
          vcMin[jz] = std::min(vcMin[jz], vcView(jnode, jz));
          vcMax[jz] = std::max(vcMax[jz], vcView(jnode, jz));
          vcAvg[jz] += vcView(jnode, jz);
        }
        ++vcSize;
      }
    }

    // Get global min/max profile
    comm_.allReduceInPlace(vcMin.begin(), vcMin.end(), eckit::mpi::min());
    comm_.allReduceInPlace(vcMax.begin(), vcMax.end(), eckit::mpi::max());
    comm_.allReduceInPlace(vcAvg.begin(), vcAvg.end(), eckit::mpi::sum());
    comm_.allReduceInPlace(vcSize, eckit::mpi::sum());

    // Normalize average
    for (size_t jz = 0; jz < nz; ++jz) {
      vcAvg[jz] /= static_cast<double>(vcSize);
    }

    // Get direction
    const bool increasing = (vcAvg[0] < vcAvg[nz-1]);

    // Get constant vertical coordinate profile
    std::vector<double> vcConst(nz);
    for (size_t jz = 0; jz < nz; ++jz) {
      // Normalized vertical index
      const double alpha = static_cast<double>(jz)/static_cast<double>(nz-1);

      // Vertical coordinate
      if (increasing) {
        vcConst[jz] = alpha*vcMax[jz] + (1.0-alpha)*vcMin[jz];
      } else {
        vcConst[jz] = alpha*vcMin[jz] + (1.0-alpha)*vcMax[jz];
      }
    }

    // Create interpolation fields
    auto rowField = outerGeometryData.functionSpace().createField<int>(
      atlas::option::name(vcField.name() + "_row") | atlas::option::levels(2*nz));
    auto colField = outerGeometryData.functionSpace().createField<int>(
      atlas::option::name(vcField.name() + "_col") | atlas::option::levels(2*nz));
    auto SField = outerGeometryData.functionSpace().createField<double>(
      atlas::option::name(vcField.name() + "_S") | atlas::option::levels(2*nz));
    auto rowView = make_view<int, 2>(rowField);
    auto colView = make_view<int, 2>(colField);
    auto SView = make_view<double, 2>(SField);
    rowView.assign(-1);
    colView.assign(-1);
    SView.assign(0.0);
    interpFset_.add(rowField);
    interpFset_.add(colField);
    interpFset_.add(SField);

    // Compute interpolation
    for (int jnode = 0; jnode < vcField.shape(0); ++jnode) {
      // Initialization
      size_t cz = 0;
      size_t iz = 0;

      for (size_t jz = 0; jz < nz; ++jz) {
        // Loop over levels
        bool found = false;

        do {
          if (increasing) {
            if (vcView(jnode, jz) < vcConst[cz+1]) {
              // Value below the current constant level
              found = true;
            }
          } else {
            if (vcView(jnode, jz) > vcConst[cz+1]) {
              // Value above the current constant level
              found = true;
            }
          }

          if (cz == nz-2) {
            // Only up to the second to last level
            found = true;
          }

          if (!found) {
            // Increase index
            ++cz;
          }
        } while (!found);

        // Check interpolation index
        ASSERT(iz < 2*nz);

        // First interpolation coefficient
        rowView(jnode, iz) = cz;
        colView(jnode, iz) = jz;
        SView(jnode, iz) = (vcConst[cz+1]-vcView(jnode, jz))/(vcConst[cz+1]-vcConst[cz]);
        ++iz;

        // Second interpolation coefficient
        rowView(jnode, iz) = cz+1;
        colView(jnode, iz) = jz;
        SView(jnode, iz) = (vcView(jnode, jz)-vcConst[cz])/(vcConst[cz+1]-vcConst[cz]);
        ++iz;
      }
    }
  }

  oops::Log::trace() << classname() << "::OrographicInterp done" << std::endl;
}

// -----------------------------------------------------------------------------

void OrographicInterp::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  for (const auto & var : interpVars_) {
    // Get vertical coordinate name
    const std::string vcName = params_.fieldsMetaData.value().getString(var.name() + ".vert_coord");

    // Get number of levels
    const size_t nz = var.getLevels();

    // Get field
    auto field = fset[var.name()];

    // Get views
    const auto rowView = make_view<int, 2>(interpFset_[vcName + "_row"]);
    const auto colView = make_view<int, 2>(interpFset_[vcName + "_col"]);
    const auto SView = make_view<double, 2>(interpFset_[vcName + "_S"]);
    auto view = make_view<double, 2>(field);

    for (int jnode = 0; jnode < field.shape(0); ++jnode) {
      // Copy input profile and apply weight
      std::vector<double> profile(nz);
      for (size_t jz = 0; jz < nz; ++jz) {
        profile[jz] = view(jnode, jz);
        view(jnode, jz) *= 1.0-params_.orographicFactor.value();
      }

      // Apply interpolation
      for (size_t iz = 0; iz < 2*nz; ++iz) {
        view(jnode, colView(jnode, iz)) += params_.orographicFactor.value()*
          profile[rowView(jnode, iz)]*SView(jnode, iz);
      }
    }
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void OrographicInterp::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  for (const auto & var : interpVars_) {
    // Get vertical coordinate name
    const std::string vcName = params_.fieldsMetaData.value().getString(var.name() + ".vert_coord");

    // Get number of levels
    const size_t nz = var.getLevels();

    // Get fields
    auto field = fset[var.name()];

    // Get field view
    auto view = make_view<double, 2>(field);

    // Get interpolation
    const auto rowView = make_view<int, 2>(interpFset_[vcName + "_row"]);
    const auto colView = make_view<int, 2>(interpFset_[vcName + "_col"]);
    const auto SView = make_view<double, 2>(interpFset_[vcName + "_S"]);

    for (int jnode = 0; jnode < field.shape(0); ++jnode) {
      // Copy input profile and apply weight
      std::vector<double> profile(nz);
      for (size_t jz = 0; jz < nz; ++jz) {
        profile[jz] = view(jnode, jz);
        view(jnode, jz) *= 1.0-params_.orographicFactor.value();
      }

      // Apply interpolation
      for (size_t iz = 0; iz < 2*nz; ++iz) {
        view(jnode, rowView(jnode, iz)) += params_.orographicFactor.value()*
          profile[colView(jnode, iz)]*SView(jnode, iz);
      }
    }
  }

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void OrographicInterp::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
