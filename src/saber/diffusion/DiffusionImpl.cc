/*
 * (C) Copyright 2023-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/diffusion/DiffusionImpl.h"

#include <algorithm>
#include <limits>

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace saber {
namespace diffusion {

// --------------------------------------------------------------------------------------
/// helper class to print stats about fields

class FieldStats {
 public:
  FieldStats(const atlas::Field & field, const eckit::mpi::Comm & comm) {
    bool skipZero = true;
    int count = 0;
    double sum = 0.0;
    min = std::numeric_limits<double>::max();
    max = std::numeric_limits<double>::min();

    const auto & view = atlas::array::make_view<double, 2>(field);
    for (atlas::idx_t i = 0; i < field.shape(0); i++) {
      for (atlas::idx_t lvl = 0; lvl < field.shape(1); lvl++) {
        double val = view(i, lvl);
        if (skipZero && val == 0.0) continue;
        min = std::min(min, val);
        max = std::max(max, val);
        sum += val;
        count++;
      }
    }

    comm.allReduceInPlace(count, eckit::mpi::sum());
    comm.allReduceInPlace(min, eckit::mpi::min());
    comm.allReduceInPlace(max, eckit::mpi::max());
    comm.allReduceInPlace(sum, eckit::mpi::sum());
    mean = count > 0 ? sum / count : 0.0;
  }

  void print(std::ostream &os) const {
    os << "min: " << min << "  mean: " << mean << "  max: " << max;
  }

 private:
  double min, max, mean;
};
std::ostream &operator<< (std::ostream &os, const FieldStats &fs) {fs.print(os); return os;}

// --------------------------------------------------------------------------------
/// Parse the vertical numerical scheme from the "vertical:" block of YAML and
/// apply it (plus the precomputed scales) to a Diffusion instance. The YAML keys
/// are optional; defaults preserve the pre-existing explicit-scheme behavior:
///    vertical:
///      method: explicit | implicit    # default "explicit"
///      iterations: 4                  # only used when method: implicit, must be even
///
/// The input length scale is interpreted as the Daley length scale of the
/// output kernel. createScales() applies the optional GC-half-width to Daley
/// length conversion (1/3.67) beforehand when "as gaussian: false" (the
/// default), so this routine can be shape-agnostic.
void configureDiffusion(oops::Diffusion & diffusion,
                        const atlas::FieldSet & scales,
                        const oops::OptionalParameter<eckit::LocalConfiguration> & vtConf) {
  std::string methodStr = "explicit";
  int implicitIterations = 4;
  if (vtConf.value() != boost::none) {
    methodStr = (*vtConf.value()).getString("method", "explicit");
    implicitIterations = (*vtConf.value()).getInt("iterations", 4);
  }

  if (methodStr == "explicit") {
    diffusion.setParameters(scales);
  } else if (methodStr == "implicit") {
    if (implicitIterations <= 0 || implicitIterations % 2 != 0)
      throw eckit::Exception("saber::Diffusion: vertical \"iterations\" must be a positive"
                             " even integer, got: " + std::to_string(implicitIterations));
    oops::Log::info() << "  vertical method: implicit (iterations: "
                      << implicitIterations
                      << ", length scale interpreted as Daley length)"
                      << std::endl;
    diffusion.setParameters(scales,
                            oops::Diffusion::VerticalMethod::Implicit,
                            implicitIterations);
  } else {
    throw eckit::Exception("saber::Diffusion: vertical \"method\" must be"
                           " \"explicit\" or \"implicit\", got: " + methodStr);
  }
}

// --------------------------------------------------------------------------------
// helper function to read / generate scales. The method of specifying scales is the
// same for horizontal and vertical scales

void createScales(const oops::GeometryData & geom,
                  std::queue<atlas::Field> & calibrateReadFields,
                  const std::string & scaleName,
                  const eckit::LocalConfiguration & conf,
                  atlas::FieldSet & scales) {
  const size_t levels = conf.getInt("levels", 1);

  // generate or read the length scales
  atlas::Field localScales;
  if (conf.has("fixed value")) {
    // generate scales from a fixed value
    localScales = geom.functionSpace().createField<double>(
      atlas::option::levels(levels) | atlas::option::name(scaleName));
    auto v_localScales = atlas::array::make_view<double, 2>(localScales);
    const double val = conf.getDouble("fixed value");
    v_localScales.assign(val);
    oops::Log::info() << "  fixed value: " << val << std::endl;

  } else if (conf.has("model file")) {
    oops::Log::info() << "  Reading from file" << std::endl;
    localScales = calibrateReadFields.front();
    calibrateReadFields.pop();
    localScales.haloExchange();
    localScales.rename(scaleName);
  } else {
    throw eckit::Exception("createScales: either \"fixed value\" or \"model file\""
                           " should be set.");
  }

  // apply optional masking.
  // scales are simply set to 0.0 where the mask is
  const std::string maskName = conf.getString("mask", "");
  oops::Log::info() << "  mask: " << (maskName == "" ? "NONE" : maskName) << std::endl;
  if (maskName != "") {
    const auto & v_mask = atlas::array::make_view<double, 2>(geom.getField(maskName));
    auto v_field = atlas::array::make_view<double, 2>(localScales);
    for (atlas::idx_t i = 0; i < localScales.shape(0); i++) {
      if (!v_mask(i, 0)) {
        for (size_t level = 0; level < levels; level++) {
          v_field(i, level) = 0.0;
        }
      }
    }
  }

  // apply optional conversion from GC half width to Daley length scale.
  // oops::Diffusion (both explicit and implicit methods) interprets the scale
  // it receives as the Daley length scale of the output kernel. Input to
  // saber::Diffusion, by default, assumes a Gaspari-Cohn compact half-width
  // instead; the 1/3.67 factor converts that to a Daley length, which is then
  // a shape-agnostic length parameter that works for both Gaussian-like
  // (explicit) and Matern-like (implicit) output kernels. The factor is an
  // empirical match rather than an exact conversion -- the kernel shape is
  // never actually GC, just tuned to a similar effective range.
  const bool isGaussian = conf.getBool("as gaussian", false);
  oops::Log::info() << "  scale type: " << (isGaussian ? "Gaussian / Daley" : "Gaspari-Cohn")
                      << std::endl;
  if (!isGaussian) {
    atlas::FieldSet fset;
    fset.add(localScales);
    util::multiplyFieldSet(fset, 1.0/3.67);
  }

  // all done
  scales.add(localScales);
}

// ----------------------------------------------------------------------------------
// helper function to apply the sqrt of normalization

void applyNormSqrt(const atlas::FunctionSpace & fs,
                   const Group & group,
                   atlas::Field & field) {
  // This code assumes
  // 1) the number of levels of hzNorm is greater than or equal to the number of
  //    levels in field. OR:
  // 2) hzNorm is 2D, in which case the same norm is applied to each level of field

  field.haloExchange();
  auto view = atlas::array::make_view<double, 2>(field);

  // horizontal normalization
  if (group.normalization.has("hzNorm")) {
    const atlas::Field & normHz = group.normalization["hzNorm"];
    ASSERT(normHz.shape(0) == field.shape(0));
    ASSERT(normHz.shape(1) == 1 || normHz.shape(1) >= field.shape(1));

    const auto & v_normHz = atlas::array::make_view<double, 2>(normHz);
    for (atlas::idx_t i = 0; i < fs.size(); i++) {
      for (atlas::idx_t lvl = 0; lvl < field.shape(1); lvl++) {
        const atlas::idx_t normLvl = std::min(normHz.shape(1)-1, lvl);
        view(i, lvl) *= v_normHz(i, normLvl);
      }
    }
  }

  // vertical normalization
  if (group.normalization.has("vtNorm") && field.shape(1) > 1) {
    const atlas::Field &vtNorm = group.normalization["vtNorm"];
    ASSERT(vtNorm.shape(0) == field.shape(0));
    ASSERT(vtNorm.shape(1) == field.shape(1));

    const auto & v_normVt = atlas::array::make_view<double, 2>(vtNorm);
    for (atlas::idx_t i = 0; i < fs.size(); i++) {
      for (atlas::idx_t lvl = 0; lvl < field.shape(1); lvl++) {
        view(i, lvl) *= v_normVt(i, lvl);
      }
    }
  }
}

// --------------------------------------------------------------------------------------

void addConf(const oops::OptionalParameter<eckit::LocalConfiguration> & conf,
             std::vector<std::pair<std::string, eckit::LocalConfiguration>> & confs) {
  const std::string MODEL_FILE_PARAM = "model file";
  if (conf.value() != boost::none) {
    const auto & confVal = (*conf.value());
    if (confVal.has(MODEL_FILE_PARAM)) {
      confs.push_back(std::pair<std::string, eckit::LocalConfiguration>(
        confVal.getString("model variable"),
        confVal.getSubConfiguration(MODEL_FILE_PARAM)));
    }
  }
}

// --------------------------------------------------------------------------------------

void randomize(const oops::GeometryData & geom,
               const std::vector<Group> & groups,
               oops::FieldSet3D & fset) {
  for (const auto & group : groups) {
    for (const auto & var : group.vars.variables()) {
      if (!fset.has(var)) continue;

      // create random field
      const size_t levels = fset[var].shape(1);
      atlas::FieldSet rand = util::createRandomFieldSet(
        geom.comm(), geom.functionSpace(),
        std::vector<size_t>{levels}, std::vector<std::string>{var});
      fset[var] = rand[var];

      // apply diffusion square root
      group.diffusion->multiplySqrtTL(rand);
    }
  }
}

// --------------------------------------------------------------------------------------

void multiply(const oops::GeometryData & geom,
              const std::vector<Group> & groups,
              oops::FieldSet3D & fset) {
  const atlas::FunctionSpace & fs = geom.functionSpace();

  // iterate through the list of groups
  for (const auto & group : groups) {
    // ----------------------------------------------------------------------------------

    // get the subset of fields, or create a common field if doing duplicated variable
    // strategy
    atlas::FieldSet fieldSubset;
    if (group.varDuplicated) {
      // "duplicated" multivariate strategy.
      // find the number of levels needed for a common field.
      atlas::idx_t levels = 1;
      for (const auto & var : group.vars.variables()) {
        levels = std::max(levels, fset[var].shape(1));
      }

      // create a new 3d field that is a summation of all the variable in this group
      atlas::Field common = geom.functionSpace().createField<double>(
        atlas::option::levels(levels) | atlas::option::name("COMMON"));
      fieldSubset.add(common);
      auto v_common = atlas::array::make_view<double, 2>(common);
      v_common.assign(0.0);
      for (const auto & var : group.vars.variables()) {
        const atlas::Field & field = fset[var];
        field.haloExchange();
        const auto & v_field = atlas::array::make_view<double, 2> (field);
        for (atlas::idx_t i = 0; i < fs.size(); i++) {
          for (atlas::idx_t lvl = 0; lvl < field.shape(1); lvl++) {
            v_common(i, lvl) += v_field(i, lvl);
          }
        }
      }
    } else {
      for (const auto & var : group.vars.variables()) {
        // TODO(Travis) make sure the variables in the active var list
        fieldSubset.add(fset[var]);
      }
    }

    // Do the diffusion for each field
    for (auto & field : fieldSubset) {
      atlas::FieldSet fset;
      if (group.vtDuplicated) {
        // merge vertical levels
        atlas::Field common = fs.createField<double>(atlas::option::levels(1));
        auto v_common = atlas::array::make_view<double, 2>(common);
        auto v_field = atlas::array::make_view<double, 2>(field);
        v_common.assign(0.0);
        for (atlas::idx_t i = 0; i < fs.size(); i++) {
          for (atlas::idx_t lvl = 0; lvl < field.shape(1); lvl++) {
            v_common(i, 0) += v_field(i, lvl);
          }
        }
        fset.add(common);

        // do diffusion
        applyNormSqrt(fs, group, common);
        group.diffusion->multiplySqrtAD(fset);
        group.diffusion->multiplySqrtTL(fset);
        applyNormSqrt(fs, group, common);

        // put vertical levels back to source
        for (atlas::idx_t i = 0; i < fs.size(); i++) {
          for (atlas::idx_t lvl = 0; lvl < field.shape(1); lvl++) {
            v_field(i, lvl) = v_common(i, 0);
          }
        }
      } else {
        fset.add(field);
        applyNormSqrt(fs, group, field);
        group.diffusion->multiplySqrtAD(fset);
        group.diffusion->multiplySqrtTL(fset);
        applyNormSqrt(fs, group, field);
      }
    }

    // if a duplicated var strategy, copy back to the source fields
    if (group.varDuplicated) {
      const auto & v_common = atlas::array::make_view<double, 2>(fieldSubset["COMMON"]);
      // copy back into the separate input fields
      for (const auto & var : group.vars.variables()) {
        atlas::Field & field = fset[var];
        auto v_field = atlas::array::make_view<double, 2> (field);
        for (atlas::idx_t i = 0; i < fs.size(); i++) {
          for (atlas::idx_t lvl = 0; lvl < field.shape(1); lvl++) {
            v_field(i, lvl) = v_common(i, lvl);
          }
        }
      }
    }

    // end of group loop
  }
}

// --------------------------------------------------------------------------------------

void filter(const std::vector<Group> & groups,
            oops::FieldSet3D & fset) {
  // iterate through the list of groups
  for (const auto & group : groups) {
    // get the subset of fields, or create a common field if doing duplicated variable
    // strategy
    atlas::FieldSet fieldSubset;
    for (const auto & var : group.vars.variables()) {
      // TODO(Travis) make sure the variables in the active var list
      fieldSubset.add(fset[var]);
    }
    // Do the diffusion for each field
    for (auto & field : fieldSubset) {
      atlas::FieldSet fset;
      fset.add(field);
      group.diffusion->multiplySqrtAD(fset);
      group.diffusion->multiplySqrtTL(fset);
    }
    // end of group loop
  }
}

// --------------------------------------------------------------------------------------

void read(const oops::GeometryData & geom,
          std::vector<Group> & groups,
          const DiffusionParameters & params) {
  oops::Log::info()
    << "\n==================================================================================\n"
    << " saber::Diffusion read\n"
    << "----------------------------------------------------------------------------------\n";

  const atlas::FunctionSpace & fs = geom.functionSpace();
  const DiffusionParameters::Read & readParams = *params.read.value();
  const std::vector<DiffusionParameters::Group> & paramGroups = readParams.groups.value();

  // process all of the listed groups
  int groupCount = 1;
  for (const DiffusionParameters::Group & groupConf : paramGroups) {
    if (groupCount > 1) {
      oops::Log::info()
        << "\n----------------------------------------------------------------------------------\n";
    }
    oops::Log::info() << "\ngroup " << groupCount++ << " of " << paramGroups.size() <<std::endl;

    Group & group = groups.emplace_back();
    group.diffusion.reset(new oops::Diffusion(geom, oops::Diffusion::calculateDerivedGeom(geom)));

    // sanity check on input config
    if (groupConf.vertical.value() == boost::none && groupConf.horizontal.value() == boost::none) {
      throw eckit::Exception("Diffusion: At least 1 of \"horizontal\" or \"vertical\" must be"
                             " specified in a read group");
    }

    // get applicable vars
    if (groupConf.variables.value() != boost::none) {
      group.vars = oops::Variables(*groupConf.variables.value());
    }
    oops::Log::info() << "variables: " << group.vars.variables() << std::endl;

    // multivariate strategy
    group.varDuplicated = groupConf.multivariateStrategy.value() ==
                          DiffusionParameters::Strategy::DUPLICATED;
    oops::Log::info() << "multivariate strategy: "
                      << (group.varDuplicated ? "duplicated" : "univariate") << std::endl;

    // read scales and normalization.
    //
    // TODO(Travis) I really should be reading in the diffusion constants and number of
    // iterations, instead of the scales. Need to change oops::Diffusion for that
    // TODO(Travis) I don't like having to specify the number of levels when doing a read.
    // modify util::readFieldSet to not need this?
    atlas::FieldSet diffusionParams;
    if (groupConf.horizontal.value() != boost::none) {
      // read in the horizontal parameters
      oops::Log::info() << "Reading horizontal parameters" << std::endl;
      const eckit::LocalConfiguration hzConf = *groupConf.horizontal.value();

      const size_t levels = hzConf.getInt("levels", 1);
      oops::Log::info() << "  levels: " << levels << std::endl;

      atlas::FieldSet hzParams;
      oops::Log::info() << "  ";
      util::readFieldSet(geom.comm(), fs,
        std::vector<size_t>{levels, levels},
        std::vector<std::string>{"hzScales", "hzNorm"},
        hzConf,
        hzParams);
      hzParams.haloExchange();
      util::shareFields(hzParams, diffusionParams);
      group.normalization.add(hzParams["hzNorm"]);

      oops::Log::info() << "  hzScales:  " << FieldStats(hzParams["hzScales"], geom.comm())
                        << std::endl;
      oops::Log::info() << "  hzNorm:    " << FieldStats(hzParams["hzNorm"], geom.comm())
                        << std::endl;
    }

    if (groupConf.vertical.value() != boost::none) {
      // read in the vertical parameters
      oops::Log::info() << "Reading vertical parameters" << std::endl;
      const eckit::LocalConfiguration vtConf = *groupConf.vertical.value();

      const size_t levels = vtConf.getInt("levels", 1);
      oops::Log::info() << "  levels: " << levels << std::endl;

      std::string verticalStrategy = vtConf.getString("strategy", "normal");
      group.vtDuplicated = (verticalStrategy == "duplicated");
      oops::Log::info() << "  vertical strategy: " << verticalStrategy << std::endl;

      if (!group.vtDuplicated) {
        atlas::FieldSet vtParams;
        oops::Log::info() << "  ";
        util::readFieldSet(geom.comm(), fs,
          std::vector<size_t>{levels, levels},
          std::vector<std::string>{"vtScales", "vtNorm"},
          *groupConf.vertical.value(),
          vtParams);
        vtParams.haloExchange();
        util::shareFields(vtParams, diffusionParams);
        group.normalization.add(diffusionParams["vtNorm"]);

        oops::Log::info() << "  vtScales:  " << FieldStats(vtParams["vtScales"], geom.comm())
                          << std::endl;
        oops::Log::info() << "  vtNorm:    " << FieldStats(vtParams["vtNorm"], geom.comm())
                          << std::endl;
      }
    }

    // set the diffusion scales (and select vertical scheme if configured)
    configureDiffusion(*group.diffusion, diffusionParams, groupConf.vertical);
  }
  oops::Log::info()
    << "==================================================================================\n\n";
}

// --------------------------------------------------------------------------------------

std::vector<std::pair<std::string, eckit::LocalConfiguration>> getReadConfs(
  const DiffusionParameters & params) {
  // Set the config for the model files that are to be read during calibration.
  std::vector<std::pair<std::string, eckit::LocalConfiguration>> confs;

  if (params.calibration.value() != boost::none) {
    const auto & calibrationConf = (*params.calibration.value());

    // for each diffusion group...
    for (const auto & group : calibrationConf.groups.value()) {
      addConf(group.horizontal, confs);
      addConf(group.vertical, confs);
    }
  }
  return confs;
}

// --------------------------------------------------------------------------------------

void setReadFields(const std::vector<oops::FieldSet3D> & fvec,
                   std::queue<atlas::Field> & calibrateReadFields) {
  // Save the atlas fields needed for calibration. Note that for simplicity we
  // cram these fields into a queue because the fieldsets in fvec are going to
  // be in the read same order as the the eckit configuration used later in
  // directCalibration(). I'm lazy (but efficient?)
  for (const auto & fset : fvec) {
    calibrateReadFields.push(fset[fset.name()]);
  }
}

// --------------------------------------------------------------------------------------

void directCalibration(const oops::GeometryData & geom,
                       std::vector<Group> & groups,
                       std::queue<atlas::Field> & calibrateReadFields,
                       const DiffusionParameters & params) {
  oops::Log::info()
    << "\n==================================================================================\n"
    << " saber::Diffusion calibration\n"
    << "----------------------------------------------------------------------------------\n";

  const atlas::FunctionSpace & fs = geom.functionSpace();

  // get relevant configuration
  const DiffusionParameters::Calibration & calibrationParams = *params.calibration.value();
  const int randomizationIterations = calibrationParams.normalizationIterations.value();
  const std::vector<DiffusionParameters::Group> & paramGroups = calibrationParams.groups.value();

  // process all of the listed groups
  int groupCount = 1;
  for (const DiffusionParameters::Group & groupConf : paramGroups) {
    if (groupCount > 1) {
      oops::Log::info()
        << "\n----------------------------------------------------------------------------------\n";
    }
    oops::Log::info() << "\ngroup " << groupCount++ << " of " << paramGroups.size() <<std::endl;

    Group & group = groups.emplace_back();

    // get applicable vars
    if (groupConf.variables.value() != boost::none) {
      group.vars = oops::Variables(*groupConf.variables.value());
    }

    // ------------------------------------------------------------------------------------
    // read in or generate the scales
    // ------------------------------------------------------------------------------------
    atlas::FieldSet scales;
    group.diffusion.reset(new oops::Diffusion(geom, oops::Diffusion::calculateDerivedGeom(geom)));
    // sanity check on input config
    if (groupConf.vertical.value() == boost::none &&
        groupConf.horizontal.value() == boost::none) {
      throw eckit::Exception("saber::Diffusion: At least 1 of \"horizontal\" or \"vertical\" "
                             "must be specified in a read group");
    }

    // --------------------------------------------------------------------------------

    // generate horizontal scales
    if (groupConf.horizontal.value() != boost::none) {
      oops::Log::info() << "Generating horizontal scales" << std::endl;
      const eckit::LocalConfiguration & hzConf = *groupConf.horizontal.value();
      createScales(geom, calibrateReadFields, "hzScales", hzConf, scales);
      oops::Log::info() << "  hzScales:  " << FieldStats(scales["hzScales"], geom.comm())
                        << std::endl;
    }

    // create vertical scales
    if (groupConf.vertical.value() != boost::none) {
      oops::Log::info() << "Generating vertical scales" << std::endl;
      const eckit::LocalConfiguration & vtConf = *groupConf.vertical.value();
      createScales(geom, calibrateReadFields, "vtScales", vtConf, scales);
      ASSERT(scales["vtScales"].shape(1) > 1);
      oops::Log::info() << "  vtScales:  " << FieldStats(scales["vtScales"], geom.comm())
                        << std::endl;
    }

    // done specifying scales. Pass them to oops::Diffusion
    // (and select vertical scheme if configured)
    configureDiffusion(*group.diffusion, scales, groupConf.vertical);

    // ------------------------------------------------------------------------------------
    // Calculate horizontal normalization. This is done using a randomization method that
    // will eventually converge to the true normalization coefficients, given enough
    // iterations. The process is to 1) create a random field 2) apply the sqrt of
    // diffusion 3) keep track of the running variance 4) repeat, and when done, the
    // normalization is a function of this variance.
    // ------------------------------------------------------------------------------------
    if (groupConf.horizontal.value() != boost::none) {
      oops::Log::info() << "Calculating horizontal normalization...\n";

      const size_t levels = scales["hzScales"].shape(1);

      auto normHz = fs.createField<double>(
        atlas::option::levels(levels) |
        atlas::option::name("hzNorm"));
      group.normalization.add(normHz);
      auto v_normHz = atlas::array::make_view<double, 2>(normHz);
      v_normHz.assign(1.0);

      // fields that are needed to keep a running variance calculation
      atlas::Field s = fs.createField<double>(atlas::option::levels(levels));
      atlas::Field m = fs.createField<double>(atlas::option::levels(levels));
      auto v_s = atlas::array::make_view<double, 2>(s);
      auto v_m = atlas::array::make_view<double, 2>(m);
      v_s.assign(0.0);
      v_m.assign(0.0);

      // Perform multiple iterations of calculating the variance of the diffusion operator
      // when random vectors are supplied
      oops::Log::info() << "  randomization iterations: " << randomizationIterations << std::endl;
      oops::Log::info() << "  status:  0%"<< std::endl;
      const int n10pct = std::floor(randomizationIterations/10.0);
      for (int itr = 1; itr <= randomizationIterations; itr++) {
        if ( itr % n10pct == 0 ) {
          oops::Log::info() << "          " << 10*itr/n10pct << "%" << std::endl;
        }

        // generate random vector
        atlas::FieldSet rand = util::createRandomFieldSet(geom.comm(), fs,
          std::vector<size_t>{levels}, std::vector<std::string>{"rand"});

        // apply sqrt of horizontal diffusion
        group.diffusion->multiplySqrtTL(rand, oops::Diffusion::Mode::HorizontalOnly);

        // keep track of the stats needed for a running variance calculation
        // (Welford 1962 algorithm)
        const auto & v_rand = atlas::array::make_view<double, 2>(rand["rand"]);
        double new_m;
        for (atlas::idx_t i = 0; i < fs.size(); i++) {
          for (size_t lvl = 0; lvl < levels; lvl++) {
            const double m = v_m(i, lvl);
            const double f = v_rand(i, lvl);
            new_m = m + (f - m) / itr;
            v_s(i, lvl) += (f - m)*(f - new_m);
            v_m(i, lvl) = new_m;
          }
        }
      }  // done with randomization iterations

      // calculate final normalization coefficients
      for (atlas::idx_t i = 0; i < fs.size(); i++) {
        for (size_t lvl = 0; lvl < levels; lvl++) {
          if (v_s(i, lvl) > 0.0) {
            v_normHz(i, lvl) = 1.0 / sqrt(v_s(i, lvl) / (randomizationIterations-1));
          }
        }
      }
    }

    // ------------------------------------------------------------------------------------
    // Calculate vertical normalization. This is done using a brute force method. A dirac
    // is created at each level, vertical diffusion is applied, and the normalization is
    // calculated from that
    // ------------------------------------------------------------------------------------
    if (groupConf.vertical.value() != boost::none) {
      oops::Log::info() << "Calculating vertical normalization...\n";
      const size_t levels = scales["vtScales"].shape(1);

      auto normVt = fs.createField<double>(
        atlas::option::levels(levels) |
        atlas::option::name("vtNorm"));
      group.normalization.add(normVt);
      auto v_normVt = atlas::array::make_view<double, 2>(normVt);
      v_normVt.assign(1.0);

      atlas::FieldSet normLvl = util::createFieldSet(fs,
        std::vector<size_t>{levels}, std::vector<std::string>{"norm"});
      auto v_normLvl = atlas::array::make_view<double, 2>(normLvl["norm"]);

      // for each level
      for (size_t lvl = 0; lvl < levels; lvl++) {
        // assign 0 everywhere except 1.0 for the current level
        v_normLvl.assign(0.0);
        for (atlas::idx_t i = 0; i < fs.size(); i++) {
          v_normLvl(i, lvl) = 1.0;
        }

        // apply vertical diffusion
        group.diffusion->multiplySqrtAD(normLvl, oops::Diffusion::Mode::VerticalOnly);
        group.diffusion->multiplySqrtTL(normLvl, oops::Diffusion::Mode::VerticalOnly);

        // save the normalization coefficients for this level
        for (atlas::idx_t i = 0; i < fs.size(); i++) {
          if (v_normLvl(i, lvl) > 0.0) {
            v_normVt(i, lvl) = 1.0 / sqrt(v_normLvl(i, lvl));
          }
        }
      }
    }

    // ------------------------------------------------------------------------------------
    // write out calibration parameters (scales and normalization)
    // TODO(travis) ideally the diffusion coefficients and number of iterations should be
    // directly saved, instead of re-calculated every time from the given scales
    // ------------------------------------------------------------------------------------
    if (groupConf.write.value() != boost::none) {
      atlas::FieldSet writeFields = util::shareFields(scales);
      util::shareFields(group.normalization, writeFields);
      util::writeFieldSet(geom.comm(), *groupConf.write.value(), writeFields);
    }
  }
  oops::Log::info()
    << "==================================================================================\n\n";
}

// --------------------------------------------------------------------------------------

}  // namespace diffusion
}  // namespace saber
