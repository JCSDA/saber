/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/generic/WriteFields.h"

#include <sys/stat.h>

#include <sstream>
#include <vector>

#include "oops/base/FieldSet3D.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<WriteFields> makerWriteFields_("write fields");

// -----------------------------------------------------------------------------

void WriteFields::writeToFile(const atlas::FieldSet & fset,
                              const std::string & description,
                              size_t & count) const {
  const auto & fieldNames = params_.fieldNames.value();

  // Select the fields to write out.
  // If the `fieldNames` list is empty, write out all fields in `fset`.
  atlas::FieldSet fsetwrite;
  for (const auto & field : fset) {
    if (fieldNames.empty() ||
        std::find(fieldNames.begin(), fieldNames.end(), field.name()) != fieldNames.end()) {
      fsetwrite.add(field);
    }
  }

  // Check there is at least one field to write out.
  if (fsetwrite.size() == 0) {
    oops::Log::warning() << "No fields present. File will not be written." << std::endl;
    return;
  }

  // Produce a filename and increment the relevant counter.
  std::stringstream filepath;
  filepath << params_.outputPath.value() << "/"
           << "fields_" << params_.label.value() << "_" << description << "_" << count++;
  const std::string filepathnc = filepath.str() + ".nc";

  eckit::LocalConfiguration conf;
  conf.set("filepath", filepath.str());

  // todo: replace the use of `stat` with `std::filesystem::exists`
  // when compilers have been upgraded.
  struct stat buffer;
  if (stat(filepathnc.c_str(), &buffer) == 0) {
    oops::Log::warning() << "File " << filepathnc << " already exists. Overwriting." << std::endl;
  }

  // Write field set to file.
  util::writeFieldSet(innerGeometryData_.comm(), conf, fsetwrite);

  // Output for validation.
  oops::Log::test() << "Wrote file " << filepathnc << " with fields:" << std::endl;
  for (const auto & field : fsetwrite) {
    oops::Log::test() << "  " << field.name() << ": "
                      << util::normField(field, innerGeometryData_.comm()) << std::endl;
  }
  oops::Log::test() << std::endl;
}

// -----------------------------------------------------------------------------

WriteFields::WriteFields(const oops::GeometryData & outerGeometryData,
                         const oops::Variables & outerVars,
                         const eckit::Configuration & covarConfig,
                         const Parameters_ & params,
                         const oops::FieldSet3D & xb,
                         const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params),
    innerGeometryData_(outerGeometryData),
    innerVars_(outerVars),
    params_(params),
    count_xb_(0),
    count_fg_(0),
    count_multiply_(0),
    count_multiplyad_(0),
    count_leftinversemultiply_(0)
{
  oops::Log::trace() << classname() << "::WriteFields starting" << std::endl;

  if (params_.writeXb) {
    writeToFile(xb.fieldSet(), "xb", count_xb_);
  }

  if (params_.writeFg) {
    writeToFile(fg.fieldSet(), "fg", count_fg_);
  }

  oops::Log::trace() << classname() << "::WriteFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void WriteFields::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  if (params_.writeMultiply) {
    writeToFile(fset, "multiply", count_multiply_);
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void WriteFields::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  if (params_.writeMultiplyAD) {
    writeToFile(fset, "multiplyad", count_multiplyad_);
  }

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void WriteFields::leftInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  if (params_.writeLeftInverseMultiply) {
    writeToFile(fset, "leftinversemultiply", count_leftinversemultiply_);
  }

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void WriteFields::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
