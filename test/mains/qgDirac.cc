/*
 * (C) Copyright 2017-2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/QgTraits.h"
#include "oops/runs/Dirac.h"
#include "oops/runs/Run.h"
#include "saber/oops/instantiateCovarFactory.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  saber::instantiateCovarFactory<qg::QgTraits>();
  oops::Dirac<qg::QgTraits> dir;
  return run.execute(dir);
}
