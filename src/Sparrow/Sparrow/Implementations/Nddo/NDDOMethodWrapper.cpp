/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "NDDOMethodWrapper.h"
#include <Core/Exceptions.h>
#include <Sparrow/Implementations/Nddo/Utils/DipoleUtils/NDDODipoleMatrixCalculator.h>
#include <Sparrow/Implementations/Nddo/Utils/DipoleUtils/NDDODipoleMomentCalculator.h>
#include <Sparrow/StatesHandling/NDDOStatesHandler.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/MethodEssentials/Methods/SCFMethod.h>
#include <Utils/MethodEssentials/util/SpinAdaptedMatrix.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingPopulator.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Sparrow {

NDDOMethodWrapper::NDDOMethodWrapper() {
  this->statesHandler_ = std::make_unique<NDDOStatesHandler>(*this);
}

NDDOMethodWrapper::~NDDOMethodWrapper() = default;

Utils::PropertyList NDDOMethodWrapper::possibleProperties() const {
  return Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Dipole |
         Utils::Property::DipoleGradient | Utils::Property::DipoleMatrixAO | Utils::Property::DipoleMatrixMO |
         Utils::Property::OneElectronMatrix | Utils::Property::TwoElectronMatrix;
}

Utils::Results NDDOMethodWrapper::assembleResults(const std::string& description) const {
  bool oneElectronMatrixRequired = requiredProperties_.containsSubSet(Utils::Property::OneElectronMatrix);
  bool twoElectronMatrixRequired = requiredProperties_.containsSubSet(Utils::Property::TwoElectronMatrix);
  bool AODipoleMatrixRequired = requiredProperties_.containsSubSet(Utils::Property::DipoleMatrixAO);
  bool MODipoleMatrixRequired = requiredProperties_.containsSubSet(Utils::Property::DipoleMatrixMO);
  auto results = GenericMethodWrapper::assembleResults(description);

  auto dipoleEvaluationCoordinate = Utils::Position::Zero();

  if (oneElectronMatrixRequired) {
    results.setOneElectronMatrix(getOneElectronMatrix());
  }
  if (twoElectronMatrixRequired) {
    results.setTwoElectronMatrix(getTwoElectronMatrix());
  }
  if (AODipoleMatrixRequired) {
    if (!dipoleMatrixCalculator_->isValid()) {
      dipoleMatrixCalculator_->fillDipoleMatrix(dipoleEvaluationCoordinate);
    }
    results.setAODipoleMatrix(dipoleMatrixCalculator_->getAODipoleMatrix());
  }
  if (MODipoleMatrixRequired) {
    if (!dipoleMatrixCalculator_->isValid()) {
      dipoleMatrixCalculator_->fillDipoleMatrix(dipoleEvaluationCoordinate);
    }
    results.setMODipoleMatrix(dipoleMatrixCalculator_->getMODipoleMatrix());
  }
  return results;
}

void NDDOMethodWrapper::applySettings(std::unique_ptr<Utils::Settings>& settings, Utils::SCFMethod& method) {
  if (settings->check()) {
    bool isUnrestricted = settings->getBool(Utils::SettingsNames::unrestrictedCalculation);
    int molecularCharge = settings->getInt(Utils::SettingsNames::molecularCharge);
    int spinMultiplicity = settings->getInt(Utils::SettingsNames::spinMultiplicity);
    double selfConsistenceCriterion = settings->getDouble(Utils::SettingsNames::selfConsistanceCriterion);
    int maxIterations = settings->getInt(Utils::SettingsNames::maxIterations);
    auto scfMixerName = settings->getString(Utils::SettingsNames::mixer);
    auto scfMixerType = Utils::UniversalSettings::SettingPopulator::stringToSCFMixer(scfMixerName);
    auto logVerbosity = settings->getString(Utils::SettingsNames::loggerVerbosity);

    method.setUnrestrictedCalculation(isUnrestricted);
    method.setMolecularCharge(molecularCharge);
    method.setSpinMultiplicity(spinMultiplicity);
    method.setConvergenceCriteria(selfConsistenceCriterion);
    method.setMaxIterations(maxIterations);
    method.setScfMixer(scfMixerType);
    method.startLogger(logVerbosity);
  }
  else {
    throw Core::InitializationException("settings invalid!");
  }
}

} // namespace Sparrow
} // namespace Scine
