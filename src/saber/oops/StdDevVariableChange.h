/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_STDDEVVARIABLECHANGE_H_
#define SABER_OOPS_STDDEVVARIABLECHANGE_H_

#include <memory>
#include <string>
#include <vector>
#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"

#include "oops/base/LinearVariableChangeBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"

#include "saber/oops/OoBump.h"
#include "saber/oops/ParametersBUMP.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

/// Derived class of generic variable transform for statistical

template <typename MODEL>
class StdDevVariableChange : public oops::LinearVariableChangeBase<MODEL> {
  typedef oops::Geometry<MODEL>  Geometry_;
  typedef oops::Increment<MODEL> Increment_;
  typedef OoBump<MODEL>          OoBump_;
  typedef oops::State<MODEL>     State_;
  typedef oops::State4D<MODEL>   State4D_;
  typedef ParametersBUMP<MODEL>  Parameters_;

 public:
  static const std::string classname() {return "saber::StdDevVariableChange";}

  StdDevVariableChange(const State_ &, const State_ &,
                      const Geometry_ &, const eckit::Configuration &);
  virtual ~StdDevVariableChange();

  void multiply(const Increment_ &, Increment_ &) const override;
  void multiplyInverse(const Increment_ &, Increment_ &) const override;
  void multiplyAD(const Increment_ &, Increment_ &) const override;
  void multiplyInverseAD(const Increment_ &, Increment_ &) const override;

 private:
  const std::vector<std::string> pathStrList;

  void print(std::ostream &) const override;
  std::string ExtractFilename(std::string const pathString);
  void ExtractModelVarForCalc(std::string const fileString,
                              std::string& modelVarToCalcString,
                              std::string& varRegrByString);
  void VariableLists(const std::vector<std::string> pathStrList,
                     std::vector<std::string>& modelVarToCalcList,
                     std::vector<std::string>& varRegrByList);

  std::unique_ptr<OoBump_> ooBump_;

  // StatsVarData populate(const varin_ , const varout_, const eckit::Configuration);
};

template<typename MODEL>
StdDevVariableChange<MODEL>::StdDevVariableChange(const State_ & xb, const State_ &,
                                                const Geometry_ & resol,
                                                const eckit::Configuration & conf)
  : oops::LinearVariableChangeBase<MODEL>(conf), ooBump_()
{
  oops::Log::trace() << "StdDevVariableChange<MODEL>::StdDevVariableChange starting" << std::endl;

// Setup variables
  const eckit::LocalConfiguration varConfig(conf, "inputVariables");
  const oops::Variables vars(varConfig);

// Setup timeslots
  std::vector<util::DateTime> timeslots;
  timeslots.push_back(xb.validTime());

// Setup parameters
  Parameters_ param(resol, vars, timeslots, conf);

// Transfer OoBump pointer
  ooBump_.reset(new OoBump_(param.getOoBump()));

  oops::Log::trace() << "StdDevVariableChange<MODEL>::StdDevVariableChange done" << std::endl;
}

template<typename MODEL>
StdDevVariableChange<MODEL>::~StdDevVariableChange() {
  oops::Log::trace() << "StdDevVariableChange<MODEL>::~StdDevVariableChange starting" << std::endl;
  util::Timer timer(classname(), "~StdDevVariableChange");
  oops::Log::trace() << "StdDevVariableChange<MODEL>::~StdDevVariableChange done" << std::endl;
}

template<typename MODEL>
void StdDevVariableChange<MODEL>::multiply(const Increment_ & in, Increment_ & out) const {
  oops::Log::trace() << "StdDevVariableChange<MODEL>::multiply starting" << std::endl;

  ooBump_->multiplyStdDev(in, out);

  oops::Log::trace() << "StdDevVariableChange<MODEL>::multiply done" << std::endl;
}

template<typename MODEL>
void StdDevVariableChange<MODEL>::multiplyInverse(const Increment_ & in, Increment_ & out) const {
  oops::Log::trace() << "StdDevVariableChange<MODEL>::multiplyInverse starting" << std::endl;

  ooBump_->multiplyStdDevInv(in, out);

  oops::Log::trace() << "StdDevVariableChange<MODEL>::multiplyInverse done" << std::endl;
}

template<typename MODEL>
void StdDevVariableChange<MODEL>::multiplyAD(const Increment_ & in, Increment_ & out) const {
  oops::Log::trace() << "StdDevVariableChange<MODEL>::multiplyAD starting" << std::endl;

  ooBump_->multiplyStdDev(in, out);

  oops::Log::trace() << "StdDevVariableChange<MODEL>::multiplyAD done" << std::endl;
}

template<typename MODEL>
void StdDevVariableChange<MODEL>::multiplyInverseAD(const Increment_ & in, Increment_ & out) const {
  oops::Log::trace() << "StdDevVariableChange<MODEL>::multiplyInverseAD starting" << std::endl;

  ooBump_->multiplyStdDevInv(in, out);

  oops::Log::trace() << "StdDevVariableChange<MODEL>::multiplyInverseAD done" << std::endl;
}

template<typename MODEL>
void StdDevVariableChange<MODEL>::print(std::ostream & os) const {
}

template<typename MODEL>
std::string StdDevVariableChange<MODEL>::ExtractFilename(std::string const pathString)
{
    int tempIndex;
    for (int iCharIndex = pathString.length(); iCharIndex > 0; --iCharIndex)
    {
        if (pathString.substr(iCharIndex - 1, 1) == "/")
        {
             return pathString.substr(iCharIndex, pathString.length() - iCharIndex);
        }
        tempIndex = iCharIndex;
    }
    std::cout << "ExtractFileName: tempIndex = " << tempIndex << std::endl;
    throw std::runtime_error("FileName not extracted from path " + pathString);
}

// extract {model_variable_to_be_calculated} string and {variable_regressed_by} string
// from fileString

template<typename MODEL>
void StdDevVariableChange<MODEL>::ExtractModelVarForCalc(std::string const fileString,
                                std::string& modelVarToCalcString,
                                std::string& varRegrByString)
{
// FileName is of the form {modelVarToCalcString}__{varRegrByString}__reg.nc
// find first "__" and the last "__" to find the appropriate strings
    int firstUnderscoreIndex, lastUnderscoreIndex;
    for (int iCharIndex = 0;  iCharIndex < fileString.length() - 1; ++iCharIndex)
    {
        firstUnderscoreIndex = iCharIndex;
        if (fileString.substr(iCharIndex, 2) ==  "__")
        {
          break;
        }
    }
    for (int iCharIndex = fileString.length(); iCharIndex > 0; --iCharIndex)
    {
        lastUnderscoreIndex = iCharIndex - 1;
        if (fileString.substr(iCharIndex-1, 2) == "__")
        {
          break;
        }
    }
    modelVarToCalcString = fileString.substr(0, firstUnderscoreIndex);
    varRegrByString = fileString.substr(firstUnderscoreIndex + 2,
                                        lastUnderscoreIndex - firstUnderscoreIndex - 2);
}

// append the modelVarToCalcList and varRegrByList
template<typename MODEL>
void StdDevVariableChange<MODEL>::VariableLists(const std::vector<std::string> pathStrList,
                   std::vector<std::string>& modelVarToCalcList,
                   std::vector<std::string>& varRegrByList)
{
  std::string modelVarToCalcString;
  std::string varRegrByString;
  std::string filename;
  for (auto pathString : pathStrList)
  {
     modelVarToCalcString.clear();
     varRegrByString.clear();
     filename = ExtractFilename(pathString);
     ExtractModelVarForCalc(filename, modelVarToCalcString, varRegrByString);
     modelVarToCalcList.push_back(modelVarToCalcString);
     varRegrByList.push_back(varRegrByString);
  }
}


}  // namespace saber

#endif  // SABER_OOPS_STDDEVVARIABLECHANGE_H_
