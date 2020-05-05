/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#include "oops/runs/Run.h"
#include "oops/generic/InterpolatorFactory.h"
#include "saber/interpolation/InterpolatorFactory.h"
#include "test/generic/InterpolationInterface.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  saber::InterpolatorFactory ifactory;
  test::InterpolationInterface tests(ifactory);
  return run.execute(tests);
}
