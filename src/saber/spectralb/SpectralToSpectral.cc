/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "saber/spectralb/SpectralToSpectral.h"

#include "eckit/exception/Exceptions.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<SpectralToSpectral> makerSpectralToSpectral_("spectral to spectral");

// -----------------------------------------------------------------------------

SpectralToSpectral::SpectralToSpectral(const oops::GeometryData & outerGeometryData,
                                       const oops::Variables & outerVars,
                                       const eckit::Configuration & covarConf,
                                       const SpectralToSpectralParameters & params,
                                       const oops::FieldSet3D & xb,
                                       const oops::FieldSet3D & fg) :
    SaberOuterBlockBase(params, xb.validTime()),
    innerFunctionSpace_(2 * params.inputTruncation.value() - 1),
    outerFunctionSpace_(outerGeometryData.functionSpace()),
    innerGeometryData_(innerFunctionSpace_,
                       outerGeometryData.fieldSet(),
                       outerGeometryData.levelsAreTopDown(),
                       outerGeometryData.comm()),
    activeVars_(params.activeVars.value().value_or(outerVars))
{
  oops::Log::trace() << classname() << "::SpectralToSpectral starting" << std::endl;
  oops::Log::trace() << classname() << "::SpectralToSpectral done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToSpectral::multiply(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  truncate_or_extend(innerFunctionSpace_, outerFunctionSpace_, fieldSet.fieldSet());
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToSpectral::multiplyAD(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  truncate_or_extend(outerFunctionSpace_, innerFunctionSpace_, fieldSet.fieldSet());
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToSpectral::leftInverseMultiply(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  truncate_or_extend(outerFunctionSpace_, innerFunctionSpace_, fieldSet.fieldSet());
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToSpectral::truncate_or_extend(const atlas::functionspace::Spectral & input_fspace,
                                            const atlas::functionspace::Spectral & output_fspace,
                                            atlas::FieldSet & fieldSet) const {
  // Decides whether fieldSet should be truncated or extented, based on function spaces.
  oops::Log::trace() << classname() << "::truncate_or_extend starting" << std::endl;
  const int input_truncation = input_fspace.truncation();
  const int output_truncation = output_fspace.truncation();
  if (input_truncation < output_truncation) {
    // Extend by zero-padding
    truncateAD(input_fspace, output_fspace, fieldSet);
  } else if (input_truncation > output_truncation) {
    // Truncate
    truncate(input_fspace, output_fspace, fieldSet);
  } else {
    // Same truncation, pass.
  }
  oops::Log::trace() << classname() << "::truncate_or_extend done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToSpectral::truncate(const atlas::functionspace::Spectral & input_functionspace,
                                  const atlas::functionspace::Spectral & output_functionspace,
                                  atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::truncate starting" << std::endl;

  // Initialize empty field set
  atlas::FieldSet outFieldSet;

  const auto zonal_wavenumbers_in = input_functionspace.zonal_wavenumbers();
  const auto zonal_wavenumbers_out = output_functionspace.zonal_wavenumbers();
  const atlas::idx_t nb_zonal_wavenumbers_in = zonal_wavenumbers_in.size();
  const atlas::idx_t nb_zonal_wavenumbers_out = zonal_wavenumbers_out.size();
  const std::size_t truncation_in = input_functionspace.truncation();
  const std::size_t truncation_out = output_functionspace.truncation();
  ASSERT(truncation_in > truncation_out);

  for (auto field : fieldSet) {
    if (!activeVars_.has(field.name())) {
      // Skip passive variable
      outFieldSet.add(field);
    } else {
      // active variable
      atlas::Field outField = output_functionspace.createField<double>(
                              atlas::option::name(field.name()) |
                              atlas::option::levels(field.shape(1)));
      auto inFieldView = atlas::array::make_view<double, 2>(field);
      auto outFieldView = atlas::array::make_view<double, 2>(outField);

      // We loop over indexes of the field to be truncated
      int jnode_in = 0, jnode_out = 0;
      for (int jm=0; jm < nb_zonal_wavenumbers_in; ++jm) {
        const atlas::idx_t m = zonal_wavenumbers_in(jm);
        for (std::size_t n = m; n <= truncation_in; ++n) {
          for (auto & part : {"real", "imaginary"}) {
            (void)part;  // unused
            if (jm < nb_zonal_wavenumbers_out && n <= truncation_out) {
                // In common zone, copy values.
                for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
                  outFieldView(jnode_out, jlevel) = inFieldView(jnode_in, jlevel);
                }
                ++jnode_in; ++jnode_out;
            } else {
                // In truncated zone
                ++jnode_in;
            }
          }
        }
      }
      outFieldSet.add(outField);
    }
  }
  fieldSet = outFieldSet;

  oops::Log::trace() << classname() << "::truncate done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToSpectral::truncateAD(const atlas::functionspace::Spectral & input_functionspace,
                                    const atlas::functionspace::Spectral & output_functionspace,
                                    atlas::FieldSet & fieldSet) const {
  // The adjoint of spectral truncation is extension by zero-padding.
  oops::Log::trace() << classname() << "::truncateAD starting" << std::endl;
  // Initialize empty field set
  atlas::FieldSet outFieldSet;

  const auto zonal_wavenumbers_in = input_functionspace.zonal_wavenumbers();
  const auto zonal_wavenumbers_out = output_functionspace.zonal_wavenumbers();
  const atlas::idx_t nb_zonal_wavenumbers_in = zonal_wavenumbers_in.size();
  const atlas::idx_t nb_zonal_wavenumbers_out = zonal_wavenumbers_out.size();
  const std::size_t truncation_in = input_functionspace.truncation();
  const std::size_t truncation_out = output_functionspace.truncation();
  ASSERT(truncation_in < truncation_out);

  for (auto field : fieldSet) {
    if (!activeVars_.has(field.name())) {
      // Skip passive variable
      outFieldSet.add(field);
    } else {
      // active variable
      atlas::Field outField = output_functionspace.createField<double>(
                              atlas::option::name(field.name()) |
                              atlas::option::levels(field.shape(1)));
      auto inFieldView = atlas::array::make_view<double, 2>(field);
      auto outFieldView = atlas::array::make_view<double, 2>(outField);

      // We loop over indexes of the field to be extended
      int jnode_in = 0, jnode_out = 0;
      for (int jm=0; jm < nb_zonal_wavenumbers_out; ++jm) {
        const atlas::idx_t m = zonal_wavenumbers_out(jm);
        for (std::size_t n = m; n <= truncation_out; ++n) {
          for (auto & part : {"real", "imaginary"}) {
            (void)part;  // unused
            if (jm < nb_zonal_wavenumbers_in && n <= truncation_in) {
                // In common zone, copy values.
                for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
                  outFieldView(jnode_out, jlevel) = inFieldView(jnode_in, jlevel);
                }
                ++jnode_in; ++jnode_out;
            } else {
                // In extension zone, pad with zeros.
                for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
                  outFieldView(jnode_out, jlevel) = 0.0;
                }
                ++jnode_out;
            }
          }
        }
      }
      outFieldSet.add(outField);
    }
  }
  fieldSet = outFieldSet;

  oops::Log::trace() << classname() << "::truncateAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToSpectral::print(std::ostream & os) const {
  os << classname() << " going from truncation " << innerFunctionSpace_.truncation()
                    << " to " << outerFunctionSpace_.truncation() << "," << std::endl;
  os << "from size " << innerFunctionSpace_.size() << " to "
                     << outerFunctionSpace_.size() << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber
