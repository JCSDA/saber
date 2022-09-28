/*
 * (C) Crown Copyright 2020-2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"

#include "oops/base/Variables.h"

namespace saber {

std::vector<std::size_t> createActiveVariableSizes(const oops::Variables & activeVars,
                                                   const oops::Variables & inputVars,
                                                   const std::vector<std::size_t> & variableSizes);

void applyNtimesNplus1SpectralScaling(const oops::Variables & inputNames,
                                      const oops::Variables & outputNames,
                                      const atlas::functionspace::Spectral & specFS,
                                      const atlas::idx_t & totalWavenumber,
                                      atlas::FieldSet & fSet);

atlas::Field allocateGaussUVField(
    const atlas::FunctionSpace & gaussFS,
    const oops::Variables & activeVariables,
    const std::vector<std::size_t> & activeVariableSizes);

atlas::FieldSet allocateSpectralVortDiv(
    const atlas::functionspace::Spectral & specfs,
    const oops::Variables & activeVariables,
    const std::vector<std::size_t> & activeVariableSizes);

atlas::FieldSet convertUVToFieldSet(const atlas::Field & uvField);

atlas::Field convertUVToFieldSetAD(const atlas::FieldSet& fset);

// creating Spectral functionspace
// for now assuming that all variableSizes are the same
atlas::functionspace::Spectral
    createSpectralFunctionSpace(const atlas::StructuredGrid & gaussGrid,
                                const std::vector<std::size_t> & variableSizes);

atlas::functionspace::StructuredColumns
    createGaussFunctionSpace(const atlas::StructuredGrid & gaussGrid);

atlas::FieldSet allocateGaussFieldset(
    const atlas::functionspace::StructuredColumns & gaussFunctionSpace,
    const oops::Variables & gaussNames,
    const std::vector<std::size_t> & variableSizes);

}  // namespace saber
