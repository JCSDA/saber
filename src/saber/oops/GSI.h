/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_GSI_H_
#define SABER_OOPS_GSI_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

class GSI_Parameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GSI_Parameters, Parameters)

 public:
  // Dummy parameter
  oops::RequiredParameter<int> dummyParameter{"dummy parameter", this};
};


// -----------------------------------------------------------------------------

template<typename MODEL> class GSI {
  typedef oops::Geometry<MODEL> Geometry_;

 public:
  // Constructors
  GSI(const Geometry_ &,
       const oops::Variables &,
       const GSI_Parameters &);

  // Copy
  explicit GSI(GSI &);

  // Destructor
  ~GSI();

  // C++ interfaces
  size_t getSize() {return keyGSI_.size();}
  int getKey(int igrid) const {return keyGSI_[igrid];}
  void clearKey() {keyGSI_.clear();}

 private:
  const Geometry_ resol_;
  const oops::Variables activeVars_;
  GSI_Parameters params_;
  std::vector<int> keyGSI_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
GSI<MODEL>::GSI(const Geometry_ & resol,
                const oops::Variables & activeVars,
                const GSI_Parameters & params)
  : resol_(resol), activeVars_(activeVars), params_(params), keyGSI_() {
  oops::Log::trace() << "GSI<MODEL>::GSI construction starting" << std::endl;

  // Print dummy parameter
  oops::Log::test() << "GSI dummy parameter: " << params.dummyParameter.value() << std::endl;

  // Push back dummy key
  keyGSI_.push_back(1);

  oops::Log::trace() << "GSI:GSI constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
GSI<MODEL>::GSI(GSI & other) : keyGSI_() {
  for (unsigned int jgrid = 0; jgrid < other.getSize(); ++jgrid) {
    keyGSI_.push_back(other.getKey(jgrid));
  }
  other.clearKey();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
GSI<MODEL>::~GSI() {
  for (unsigned int jgrid = 0; jgrid < keyGSI_.size(); ++jgrid) {
    // Deallocate this instance
  }
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_GSI_H_
