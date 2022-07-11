/*
 * (C) Copyright 2021-2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/QgTraits.h"
#include "oops/runs/Run.h"
#include "saber/oops/ErrorCovarianceTraining.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  saber::ErrorCovarianceTraining<qg::QgTraits> dir;
  return run.execute(dir);
}
