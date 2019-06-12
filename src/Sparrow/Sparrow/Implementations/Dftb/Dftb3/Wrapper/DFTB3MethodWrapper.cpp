/**
 * @file DFTB3MethodWrapper.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Internal Includes */
#include "DFTB3MethodWrapper.h"
#include "DFTB3Settings.h"
#include <Sparrow/Implementations/Dftb/Utils/DipoleUtils/DFTBDipoleMatrixCalculator.h>
#include <Sparrow/Implementations/Dftb/Utils/DipoleUtils/DFTBDipoleMomentCalculator.h>
/* External Includes */
#include <Core/Exceptions.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/MethodEssentials/MethodFactories/MixerFactory.h>
#include <Utils/MethodEssentials/util/MethodExceptions.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <memory>

namespace Scine {
namespace Sparrow {

DFTB3MethodWrapper::DFTB3MethodWrapper() {
  this->settings_ = std::make_unique<DFTB3Settings>();
  requiredProperties_ = Utils::Property::Energy;
  dipoleCalculator_ = std::make_unique<DFTBDipoleMomentCalculator<dftb::DFTB3>>(method_);
  dipoleMatrixCalculator_ = DFTBDipoleMatrixCalculator<dftb::DFTB3>::create(method_);
  applySettings();
};

DFTB3MethodWrapper::DFTB3MethodWrapper(const DFTB3MethodWrapper& rhs) : DFTB3MethodWrapper() {
  copyInto(*this, rhs);
}

DFTB3MethodWrapper& DFTB3MethodWrapper::operator=(const DFTB3MethodWrapper& rhs) {
  copyInto(*this, rhs);
  return *this;
}

DFTB3MethodWrapper::~DFTB3MethodWrapper() = default;

void DFTB3MethodWrapper::applySettings() {
  if (settings_->check()) {
    bool isUnrestricted = settings_->getBool(Utils::SettingsNames::unrestrictedCalculation);
    int molecularCharge = settings_->getInt(Utils::SettingsNames::molecularCharge);
    int spinMultiplicity = settings_->getInt(Utils::SettingsNames::spinMultiplicity);
    double selfConsistenceCriterion = settings_->getDouble(Utils::SettingsNames::selfConsistanceCriterion);
    int maxIterations = settings_->getInt(Utils::SettingsNames::maxIterations);
    auto scfMixerName = settings_->getString(Utils::SettingsNames::mixer);
    auto scfMixerType = Utils::UniversalSettings::SettingPopulator::stringToSCFMixer(scfMixerName);
    auto logVerbosity = settings_->getString(Utils::SettingsNames::loggerVerbosity);

    method_.setUnrestrictedCalculation(isUnrestricted);
    method_.setMolecularCharge(molecularCharge);
    method_.setSpinMultiplicity(spinMultiplicity);
    method_.setConvergenceCriteria(selfConsistenceCriterion);
    method_.setMaxIterations(maxIterations);
    method_.setScfMixer(scfMixerType);
    method_.startLogger(logVerbosity);
  }
  else {
    throw Core::InitializationException("settings invalid!");
  }
};

std::string DFTB3MethodWrapper::name() const {
  return "DFTB3";
}

void DFTB3MethodWrapper::initialize() {
  try {
    auto parameterFile = settings_->getString(Utils::SettingsNames::parameterFile);
    auto resourceDirectory = settings_->getString(Utils::SettingsNames::parameterRootDirectory);
    auto fullPathToParameters = Utils::NativeFilenames::combinePathSegments(resourceDirectory, parameterFile);
    method_.initializeFromParameterPath(fullPathToParameters);
  }
  catch (Utils::Methods::InitializationException& e) {
    throw Core::InitializationException(e.what());
  }
}

Utils::LCAOMethod& DFTB3MethodWrapper::getLCAOMethod() {
  return method_;
}

const Utils::LCAOMethod& DFTB3MethodWrapper::getLCAOMethod() const {
  return method_;
}

void DFTB3MethodWrapper::calculateImpl(Utils::derivativeType requiredDerivative) {
  method_.calculate(requiredDerivative);
}

Utils::DensityMatrix DFTB3MethodWrapper::getDensityMatrixGuess() const {
  return method_.getDensityMatrixGuess();
}

} /* namespace Sparrow */
} /* namespace Scine */
