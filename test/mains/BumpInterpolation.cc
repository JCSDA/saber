/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#include "oops/generic/InterpolatorFactory.h"
#include "oops/runs/Run.h"
#include "saber/interpolation/InterpolatorFactory.h"
#include "test/generic/InterpolationInterface.h"

#include<iostream>

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  auto ifactory = oops::interpolation::getFactoryInstance<saber::InterpolatorFactory>().create();
  std::cout << "MSM factory created " << ifactory.classname() << std::endl;
  test::InterpolationInterface tests(ifactory);
  return run.execute(tests);
}
