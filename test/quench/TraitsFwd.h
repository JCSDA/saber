/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QUENCH_TRAITSFWD_H_
#define QUENCH_TRAITSFWD_H_

#include <string>

namespace quench {

class Geometry;
class Increment;
class State;

struct Traits {
  static std::string name() {return "quench";}

  typedef quench::Geometry         Geometry;
  typedef quench::Increment        Increment;
  typedef quench::State            State;
};

}  // namespace quench

#endif  // QUENCH_TRAITSFWD_H_
