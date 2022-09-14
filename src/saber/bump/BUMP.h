/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_BUMP_BUMP_H_
#define SABER_BUMP_BUMP_H_

#include <algorithm>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

#include "saber/bump/BUMP_Parameters.h"
#include "saber/bump/type_bump.h"


namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

class BUMP {
 public:
  // Constructor
  BUMP(const eckit::mpi::Comm &,
       const atlas::FunctionSpace &,
       const atlas::FieldSet &,
       const std::vector<size_t> &,
       const oops::Variables &,
       const BUMP_Parameters &,
       const std::vector<atlas::FieldSet> &,
       const atlas::FunctionSpace & functionSpace = NULL,
       const atlas::FieldSet & fieldSet = NULL,
       const std::vector<atlas::FieldSet> & fsetVec2 = {},
       const size_t & ens1_ne_in = 0,
       const size_t & ens2_ne_in = 0);

  // Copy-constructor
  explicit BUMP(BUMP &);

  // Destructor
  ~BUMP();

  // Accessors
  const std::vector<eckit::LocalConfiguration> memberConfig1() const {return membersConfig1_;}
  const std::vector<eckit::LocalConfiguration> memberConfig2() const {return membersConfig2_;}

  // Fortran interfaces
  void addMember(const atlas::FieldSet &, const int &, const int &) const;
  void updateVbalCov(const atlas::FieldSet &, const int &) const;
  void updateVar(const atlas::FieldSet &, const int &) const;
  void updateMom(const atlas::FieldSet &, const int &, const int &) const;
  void runDrivers() const;
  void multiplyVbal(atlas::FieldSet &) const;
  void inverseMultiplyVbal(atlas::FieldSet &) const;
  void multiplyVbalAd(atlas::FieldSet &) const;
  void inverseMultiplyVbalAd(atlas::FieldSet &) const;
  void multiplyStdDev(atlas::FieldSet &) const;
  void inverseMultiplyStdDev(atlas::FieldSet &) const;
  void randomizeNicas(atlas::FieldSet &) const;
  void multiplyNicas(atlas::FieldSet &) const;
  void multiplyPsiChiToUV(atlas::FieldSet &) const;
  void multiplyPsiChiToUVAd(atlas::FieldSet &) const;
  void getParameter(const std::string &, const int &, const int &, atlas::FieldSet &) const;
  void setNcmp(const int &, const int &) const;
  void setParameter(const std::string &, const int &, const atlas::FieldSet &) const;
  void partialDealloc() const;

 private:
  BUMP_Parameters params_;
  const oops::Variables activeVars_;
  std::vector<int> keyBUMP_;
  std::vector<eckit::LocalConfiguration> membersConfig1_;
  std::vector<eckit::LocalConfiguration> membersConfig2_;
  std::vector<oops::Variables> activeVarsPerGrid_;
};

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_BUMP_BUMP_H_
