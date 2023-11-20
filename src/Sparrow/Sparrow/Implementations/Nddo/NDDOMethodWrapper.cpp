/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "NDDOMethodWrapper.h"
#include <Core/Exceptions.h>
#include <Sparrow/Implementations/Nddo/TimeDependent/LinearResponse/CISData.h>
#include <Sparrow/Implementations/Nddo/Utils/DipoleUtils/NDDODipoleMatrixCalculator.h>
#include <Sparrow/Implementations/Nddo/Utils/DipoleUtils/NDDODipoleMomentCalculator.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <Utils/Scf/MethodInterfaces/ScfMethod.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingPopulator.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Sparrow {

NDDOMethodWrapper::NDDOMethodWrapper() = default;

NDDOMethodWrapper::~NDDOMethodWrapper() = default;

Utils::PropertyList NDDOMethodWrapper::possibleProperties() const {
  auto possibleProperties = GenericMethodWrapper::possibleProperties();
  possibleProperties.addProperty(Utils::Property::DipoleGradient);
  possibleProperties.addProperty(Utils::Property::DipoleMatrixAO);
  possibleProperties.addProperty(Utils::Property::DipoleMatrixMO);
  possibleProperties.addProperty(Utils::Property::OneElectronMatrix);
  possibleProperties.addProperty(Utils::Property::TwoElectronMatrix);
  return possibleProperties;
}

void NDDOMethodWrapper::assembleResults(const std::string& description) {
  bool oneElectronMatrixRequired = requiredProperties_.containsSubSet(Utils::Property::OneElectronMatrix);
  bool twoElectronMatrixRequired = requiredProperties_.containsSubSet(Utils::Property::TwoElectronMatrix);
  bool AODipoleMatrixRequired = requiredProperties_.containsSubSet(Utils::Property::DipoleMatrixAO);
  bool MODipoleMatrixRequired = requiredProperties_.containsSubSet(Utils::Property::DipoleMatrixMO);
  GenericMethodWrapper::assembleResults(description);

  auto dipoleEvaluationCoordinate = Utils::Position::Zero();

  if (oneElectronMatrixRequired) {
    results_.set<Utils::Property::OneElectronMatrix>(getOneElectronMatrix());
  }
  if (twoElectronMatrixRequired) {
    results_.set<Utils::Property::TwoElectronMatrix>(getTwoElectronMatrix());
  }

  if (AODipoleMatrixRequired) {
    if (!dipoleMatrixCalculator_->isValid()) {
      dipoleMatrixCalculator_->fillDipoleMatrix(dipoleEvaluationCoordinate);
    }
    results_.set<Utils::Property::DipoleMatrixAO>(dipoleMatrixCalculator_->getAODipoleMatrix());
  }
  if (MODipoleMatrixRequired) {
    if (!dipoleMatrixCalculator_->isValid()) {
      dipoleMatrixCalculator_->fillDipoleMatrix(dipoleEvaluationCoordinate);
    }
    results_.set<Utils::Property::DipoleMatrixMO>(dipoleMatrixCalculator_->getMODipoleMatrix());
  }
}

void NDDOMethodWrapper::applySettings(std::unique_ptr<Utils::Settings>& settings, Utils::ScfMethod& method) {
  GenericMethodWrapper::applySettings();
  if (settings->valid()) {
    auto spinMode = Utils::SpinModeInterpreter::getSpinModeFromString(settings->getString(Utils::SettingsNames::spinMode));
    int molecularCharge = settings->getInt(Utils::SettingsNames::molecularCharge);
    int spinMultiplicity = settings->getInt(Utils::SettingsNames::spinMultiplicity);
    double selfConsistenceCriterion = settings->getDouble(Utils::SettingsNames::selfConsistenceCriterion);
    double densityRmsdThreshold = settings->getDouble(Utils::SettingsNames::densityRmsdCriterion);
    int maxScfIterations = settings->getInt(Utils::SettingsNames::maxScfIterations);
    auto scfMixerName = settings->getString(Utils::SettingsNames::mixer);
    auto scfMixerType = Utils::UniversalSettings::SettingPopulator::stringToScfMixer(scfMixerName);

    if (spinMode == Utils::SpinMode::Any) {
      if (spinMultiplicity == 1) {
        method.setUnrestrictedCalculation(false);
      }
      else {
        method.setUnrestrictedCalculation(true);
      }
    }
    if (spinMode == Utils::SpinMode::Restricted) {
      method.setUnrestrictedCalculation(false);
    }
    if (spinMode == Utils::SpinMode::Unrestricted) {
      method.setUnrestrictedCalculation(true);
    }
    method.setMolecularCharge(molecularCharge);
    method.setSpinMultiplicity(spinMultiplicity);
    method.setConvergenceCriteria({selfConsistenceCriterion, densityRmsdThreshold});
    method.setMaxIterations(maxScfIterations);
    method.setScfMixer(scfMixerType);
  }
  else {
    settings->throwIncorrectSettings();
  }
}

CISData NDDOMethodWrapper::getCISData() const {
  return getCISDataImpl();
}

bool NDDOMethodWrapper::getZPVEInclusion() const {
  return true;
}

} // namespace Sparrow
} // namespace Scine
