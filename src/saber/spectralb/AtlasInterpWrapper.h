/*
 * (C) Crown Copyright 2021-2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_SPECTRALB_ATLASINTERPWRAPPER_H_
#define SABER_SPECTRALB_ATLASINTERPWRAPPER_H_

#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "atlas/functionspace.h"
#include "atlas/grid/detail/partitioner/MatchingMeshPartitioner.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/interpolation.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/redistribution/Redistribution.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "saber/spectralb/spectralbParameters.h"

using atlas::grid::detail::partitioner::TransPartitioner;

namespace atlas {
class Field;
class FieldSet;
}

namespace saber {
namespace spectralb {

// I will need to remove the dependency util::Printable and Object Counter
// Note 1: this will assume that the output fields are initially collocated
//         have the same number of levels and have the same partitioning.
//
template<typename MODEL>
class AtlasInterpWrapper {
  typedef spectralbParameters<MODEL> Parameters_;

 public:
  static const std::string classname() {return "saber::AtlasInterpWrapper<MODEL>";}

  AtlasInterpWrapper<MODEL>(const Parameters_ &,
                            std::shared_ptr<const atlas::FieldSet> & fieldset);

  AtlasInterpWrapper<MODEL>() {}

  void execute(const atlas::Field & srcField,
               atlas::Field & targetField) const;

  void executeAdjoint(atlas::Field & srcField,
                      const atlas::Field & targetField) const;

  void execute(const atlas::FieldSet & srcFieldSet,
               atlas::FieldSet & targetFieldSet) const;

  void executeAdjoint(atlas::FieldSet & srcFieldSet,
                      const atlas::FieldSet & targetFieldSet) const;

 private:
  void print(std::ostream &) const {}
  std::vector<std::string> fieldsetNames_;
  atlas::StructuredGrid sourceGrid_;
  atlas::grid::Distribution sourceDistribution_;
  atlas::Mesh sourceMesh_;
  atlas::functionspace::StructuredColumns sourceFunctionSpace_;

  atlas::Grid outputGrid_;
  atlas::grid::MatchingMeshPartitioner targetPartitioner_;
  atlas::Mesh targetMesh_;
  atlas::FunctionSpace targetFunctionSpace_;
  // as outputFunctionSpace but with different distribution.

  atlas::FunctionSpace outputFunctionSpace_;
  atlas::Interpolation interp_;
  atlas::Redistribution redistr_;
};

}  // namespace spectralb
}  // namespace saber
//  end of Header declaration

// Header implementation
// Because this class is just a template, the method definitions must be included for the compiler.
namespace saber {
namespace spectralb {
namespace detail {

// -----------------------------------------------------------------------------
atlas::StructuredGrid createGaussGrid(const std::string & gaussGridUid) {
  // Atlas uses "F" (followed by a number # denoting wavenumber) to denote a classic Gaussian grid
  // with # wavenumbers.
  atlas::StructuredGrid gaussGrid(gaussGridUid);
  return gaussGrid;
}

atlas::grid::Distribution createSourceDistribution(const atlas::StructuredGrid & sourceGrid) {
  return atlas::grid::Distribution(sourceGrid,
                                   atlas::grid::Partitioner(new TransPartitioner()));
}

atlas::Mesh createSourceMesh(const atlas::StructuredGrid & sourceGrid,
                             const atlas::grid::Distribution & sourceDistribution) {
  return atlas::MeshGenerator("structured").generate(sourceGrid,
                                                     sourceDistribution);
}

atlas::functionspace::StructuredColumns createSourceFunctionSpace(
    const atlas::StructuredGrid & sourceGrid,
    const atlas::grid::Distribution & sourceDistribution) {
  return atlas::functionspace::StructuredColumns(sourceGrid, sourceDistribution,
                                                 atlas::util::Config("halo", 1));
}

template<typename MODEL>
atlas::Grid createOutputGrid(const spectralbParameters<MODEL> & params) {
  std::string gridName((params.outputGridUid.value() != boost::none ?
                        params.outputGridUid.value().get() :
                        params.gaussGridUid));
  return atlas::Grid(gridName);
}

atlas::grid::MatchingMeshPartitioner createTargetPartitioner(const atlas::Mesh & sourceMesh) {
  return atlas::grid::MatchingMeshPartitioner(sourceMesh);
}

atlas::Mesh createTargetMesh(const atlas::Grid & outputGrid,
                             atlas::grid::MatchingMeshPartitioner & targetPartitioner) {
  std::string meshNameType("");
  if (outputGrid->type() == "cubedsphere") {
    meshNameType = "cubedsphere_dual";
  } else {
    meshNameType = "structured";
  }
  return atlas::MeshGenerator(meshNameType).generate(outputGrid,
                                                     targetPartitioner);
}

atlas::FunctionSpace createTargetFunctionSpace(
    const atlas::Grid & outputGrid,
    atlas::grid::MatchingMeshPartitioner & targetPartitioner,
    const atlas::Mesh & targetMesh) {
  if (atlas::StructuredGrid(targetMesh.grid())) {
    return atlas::functionspace::StructuredColumns(outputGrid, targetPartitioner,
                                                   atlas::option::halo(1));
  } else if (outputGrid.name().substr(0, 2).compare("CS") == 0) {
    // assume cubed-sphere node-columns
    return atlas::functionspace::CubedSphereNodeColumns(targetMesh);
  } else {
    ABORT("AtlasInterpWrapper:Target Function space option not supported");
    return atlas::FunctionSpace();
  }
}

atlas::FunctionSpace createOutputFunctionSpace(const atlas::FieldSet & fset) {
  return fset[0].functionspace();
}

atlas::Interpolation createAtlasInterpolation(
    const atlas::functionspace::StructuredColumns & inputFS,
    const atlas::FunctionSpace & matchingFS) {

  atlas::util::Config interpConfig;
  interpConfig.set("type", "structured-linear2D");
  interpConfig.set("adjoint", "true");

  return atlas::Interpolation(interpConfig, inputFS, matchingFS);
}

atlas::Redistribution createAtlasRedistribution(
    const atlas::FunctionSpace & matchingFS,
    const atlas::FunctionSpace & outputFS) {

  return atlas::Redistribution(matchingFS, outputFS);
}


}  // namespace detail
}  // namespace spectralb
}  // namespace saber

// ------------------------------------------------------------------------------------------------
namespace saber {
namespace spectralb {
// ------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// ATLAS INTERPOLATION WRAPPER
//-------------------------------------------------------------------------------------------------
template<typename MODEL>
AtlasInterpWrapper<MODEL>::AtlasInterpWrapper(const Parameters_ & params,
                                              std::shared_ptr<const atlas::FieldSet> & fieldset) :
  fieldsetNames_((*fieldset).field_names()),
  sourceGrid_(detail::createGaussGrid(params.gaussGridUid)),
  sourceDistribution_(detail::createSourceDistribution(sourceGrid_)),
  sourceMesh_(detail::createSourceMesh(sourceGrid_, sourceDistribution_)),
  sourceFunctionSpace_(detail::createSourceFunctionSpace(sourceGrid_, sourceDistribution_)),
  outputGrid_(detail::createOutputGrid(params)),
  targetPartitioner_(detail::createTargetPartitioner(sourceMesh_)),
  targetMesh_(detail::createTargetMesh(outputGrid_, targetPartitioner_)),
  targetFunctionSpace_(
    detail::createTargetFunctionSpace(outputGrid_, targetPartitioner_, targetMesh_)),
  outputFunctionSpace_(detail::createOutputFunctionSpace(*fieldset)),
  interp_(detail::createAtlasInterpolation(sourceFunctionSpace_, targetFunctionSpace_)),
  redistr_(detail::createAtlasRedistribution(targetFunctionSpace_, outputFunctionSpace_))
{
};

//-------------------------------------------------------------------------------------------------
template<typename MODEL>
void AtlasInterpWrapper<MODEL>::execute(const atlas::Field & srcField,
                                        atlas::Field & outputField) const {
  atlas::Field srcTmp = srcField.functionspace().createField<double>
      (atlas::option::name(srcField.name()) |
       atlas::option::levels(srcField.levels()));

  auto src_v = atlas::array::make_view<double, 2>(srcField);
  auto srcTmp_v = atlas::array::make_view<double, 2>(srcTmp);
  for (atlas::idx_t t = 0; t < srcField.shape(0); ++t) {
    for (atlas::idx_t k = 0; k < srcField.shape(1); ++k) {
      srcTmp_v(t, k) = src_v(t, k);
    }
  }
  srcTmp.haloExchange();

  auto tmpField =
    targetFunctionSpace_.createField<double>(atlas::option::name(outputField.name()) |
                                             atlas::option::levels(srcField.levels()));

  auto tmp_v = atlas::array::make_view<double, 2>(tmpField);
  tmp_v.assign(0.0);

  interp_.execute(srcTmp, tmpField);

  redistr_.execute(tmpField, outputField);

  outputField.haloExchange();
}

template<typename MODEL>
void AtlasInterpWrapper<MODEL>::executeAdjoint(atlas::Field & srcField,
                                         const atlas::Field & outputField) const {
  oops::Log::trace() << "AtlasInterpWrapper<MODEL>::executeAdjoint start "
                        "srcFieldName outputFieldName" <<
                        srcField.name() << " " << outputField.name() << std::endl;

  // create a field that is a copy of target
  // take halo adjoint of temp field
  atlas::Field tmpField =
    outputField.functionspace().createField<double>(atlas::option::name(outputField.name()) |
                                                    atlas::option::levels(outputField.levels()));

  auto tmp_v = atlas::array::make_view<double, 2>(tmpField);
  auto tgt_v = atlas::array::make_view<double, 2>(outputField);
  for (atlas::idx_t t = 0; t < outputField.shape(0); ++t) {
    for (atlas::idx_t k = 0; k < outputField.shape(1); ++k) {
      tmp_v(t, k) = tgt_v(t, k);
    }
  }
  tmpField.adjointHaloExchange();

  auto tmp2Field = targetFunctionSpace_.createField<double>(
                          atlas::option::name(outputField.name()) |
                          atlas::option::levels(outputField.levels()));
  auto tmp2_v = atlas::array::make_view<double, 2>(tmp2Field);
  tmp2_v.assign(0.0);

  redistr_.execute(tmpField, tmp2Field);

  auto src_v = atlas::array::make_view<double, 2>(srcField);
  src_v.assign(0.0);

  interp_.execute_adjoint(srcField, tmp2Field);

  oops::Log::trace() << "AtlasInterpWrapper<MODEL>::executeAdjoint done" << std::endl;
}

template<typename MODEL>
void AtlasInterpWrapper<MODEL>::execute(const atlas::FieldSet & srcFieldSet,
                                        atlas::FieldSet & targetFieldSet) const {
  for (auto & srcField : srcFieldSet) {
    execute(srcField, targetFieldSet[srcField.name()]);
  }
}

template<typename MODEL>
void AtlasInterpWrapper<MODEL>::executeAdjoint(atlas::FieldSet & srcFieldSet,
                                               const atlas::FieldSet & targetFieldSet) const {
  for (auto & srcField : srcFieldSet) {
    executeAdjoint(srcField, targetFieldSet[srcField.name()]);
  }
}

}  // namespace spectralb
}  // namespace saber
// ------------------------------------------------------------------------------------------------

#endif  // SABER_SPECTRALB_ATLASINTERPWRAPPER_H_
