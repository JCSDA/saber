/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/util/Calibration.h"

#include <algorithm>
#include <vector>

#include "oops/util/AtlasArrayUtil.h"

namespace util {

void createCalibrationNetCDFHeaderInput(const eckit::LocalConfiguration & conf,
                                        const std::string & statsType,
                                        const std::string & binType,
                                        const oops::Variables & vars,
                                        const bool & doingCalibration,
                                        eckit::LocalConfiguration & netCDFConf) {
  if (doingCalibration) {
    int fldindx(0);
    for (const auto & var : vars) {
      const std::string varName = var.name();
      const std::int32_t modellevels = var.getLevels();
      // replacing spaces for underscores in name
      std::string netCDFShortName = statsType;
      std::replace(netCDFShortName.begin(), netCDFShortName.end(), ' ', '_');
      netCDFShortName.append("_" + std::to_string(fldindx));
      const std::string longname = statsType + " of " + varName + " and " + varName;

      util::setAttribute<std::string>(netCDFConf, netCDFShortName, "long_name",
                                      "string", longname);
      util::setAttribute<std::string>(netCDFConf, netCDFShortName, "statistics_type",
                                      "string", statsType);
      util::setAttribute<std::string>(netCDFConf, netCDFShortName, "binning_type",
                                      "string", binType);
      util::setAttribute<std::string>(netCDFConf, netCDFShortName, "variable_name_1",
                                      "string", varName);
      util::setAttribute<std::string>(netCDFConf, netCDFShortName, "variable_name_2",
                                      "string", varName);
      util::setAttribute<std::int32_t>(netCDFConf, netCDFShortName, "levels_1",
                                      "int32", modellevels);
      util::setAttribute<std::int32_t>(netCDFConf, netCDFShortName, "levels_2",
                                      "int32", modellevels);
      ++fldindx;
    }

    if (conf.has("additional cross covariances")) {
      std::vector<eckit::LocalConfiguration> aconfs =
        conf.getSubConfigurations("additional cross covariances");

      for (auto & conf : aconfs) {
        const std::string var1 = conf.getString("variable 1");
        const std::string var2 = conf.getString("variable 2");
        // check that var1 and var2 are included in innerVars
        const std::int32_t levels1 = vars[var1].getLevels();
        const std::int32_t levels2 = vars[var2].getLevels();
        std::string netCDFShortName = statsType;
        std::replace(netCDFShortName.begin(), netCDFShortName.end(), ' ', '_');
        netCDFShortName.append("_" + std::to_string(fldindx));

        const std::string longname = statsType + " of " + var1 + " and " + var2;

        util::setAttribute<std::string>(netCDFConf, netCDFShortName, "long_name",
                                        "string", longname);
        util::setAttribute<std::string>(netCDFConf, netCDFShortName, "statistics_type",
                                        "string", statsType);
        util::setAttribute<std::string>(netCDFConf, netCDFShortName, "binning_type",
                                        "string", binType);
        util::setAttribute<std::string>(netCDFConf, netCDFShortName, "variable_name_1",
                                        "string", var1);
        util::setAttribute<std::string>(netCDFConf, netCDFShortName, "variable_name_2",
                                        "string", var2);
        util::setAttribute<std::int32_t>(netCDFConf, netCDFShortName, "levels_1",
                                        "int32", levels1);
        util::setAttribute<std::int32_t>(netCDFConf, netCDFShortName, "levels_2",
                                        "int32", levels2);
        ++fldindx;
      }
    }
  }
}


}  // namespace util
