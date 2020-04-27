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
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
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
  if (settings->check()) {
    bool isUnrestricted = settings->getBool(Utils::SettingsNames::unrestrictedCalculation);
    int molecularCharge = settings->getInt(Utils::SettingsNames::molecularCharge);
    int spinMultiplicity = settings->getInt(Utils::SettingsNames::spinMultiplicity);
    double selfConsistenceCriterion = settings->getDouble(Utils::SettingsNames::selfConsistanceCriterion);
    int maxIterations = settings->getInt(Utils::SettingsNames::maxIterations);
    auto scfMixerName = settings->getString(Utils::SettingsNames::mixer);
    auto scfMixerType = Utils::UniversalSettings::SettingPopulator::stringToScfMixer(scfMixerName);

    method.setUnrestrictedCalculation(isUnrestricted);
    method.setMolecularCharge(molecularCharge);
    method.setSpinMultiplicity(spinMultiplicity);
    method.setConvergenceCriteria(selfConsistenceCriterion);
    method.setMaxIterations(maxIterations);
    method.setScfMixer(scfMixerType);

    // After call to verifyPesValidity, update settings if they were internally changed due to
    // charge/spin multiplicity mismatch if not empty
    if (getLcaoMethod().getNumberAtomicOrbitals() != 0) {
      method.verifyPesValidity();
      settings->modifyBool(Utils::SettingsNames::unrestrictedCalculation, method.unrestrictedCalculationRunning());
      settings->modifyInt(Utils::SettingsNames::molecularCharge, method.getMolecularCharge());
      settings->modifyInt(Utils::SettingsNames::spinMultiplicity, method.spinMultiplicity());
    }
  }
  else {
    throw Core::InitializationException("settings invalid!");
  }
}

bool NDDOMethodWrapper::getZPVEInclusion() const {
  return true;
}

} // namespace Sparrow
} // namespace Scine
