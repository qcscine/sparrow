/**
 * @file DFTB2MethodWrapper.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Internal Includes */
#include "DFTB2MethodWrapper.h"
#include "DFTB2Settings.h"
#include <Sparrow/Implementations/Dftb/Utils/DipoleUtils/DFTBDipoleMatrixCalculator.h>
#include <Sparrow/Implementations/Dftb/Utils/DipoleUtils/DFTBDipoleMomentCalculator.h>
/* External Includes */
#include <Core/Exceptions.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/Scf/MethodExceptions.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <memory>

namespace Scine {
namespace Sparrow {

DFTB2MethodWrapper::DFTB2MethodWrapper() {
  this->settings_ = std::make_unique<DFTB2Settings>();
  requiredProperties_ = Utils::Property::Energy;
  dipoleCalculator_ = std::make_unique<DFTBDipoleMomentCalculator<dftb::DFTB2>>(method_);
  dipoleMatrixCalculator_ = DFTBDipoleMatrixCalculator<dftb::DFTB2>::create(method_);
  applySettings();
}

DFTB2MethodWrapper::DFTB2MethodWrapper(const DFTB2MethodWrapper& rhs) : DFTB2MethodWrapper() {
  copyInto(*this, rhs);
}

DFTB2MethodWrapper& DFTB2MethodWrapper::operator=(const DFTB2MethodWrapper& rhs) {
  copyInto(*this, rhs);
  return *this;
}

DFTB2MethodWrapper::~DFTB2MethodWrapper() = default;

void DFTB2MethodWrapper::applySettings() {
  GenericMethodWrapper::applySettings();
  if (settings_->check()) {
    bool isUnrestricted = settings_->getBool(Utils::SettingsNames::unrestrictedCalculation);
    int molecularCharge = settings_->getInt(Utils::SettingsNames::molecularCharge);
    int spinMultiplicity = settings_->getInt(Utils::SettingsNames::spinMultiplicity);
    double selfConsistenceCriterion = settings_->getDouble(Utils::SettingsNames::selfConsistanceCriterion);
    int maxIterations = settings_->getInt(Utils::SettingsNames::maxIterations);
    auto scfMixerName = settings_->getString(Utils::SettingsNames::mixer);
    auto scfMixerType = Utils::UniversalSettings::SettingPopulator::stringToScfMixer(std::move(scfMixerName));

    method_.setUnrestrictedCalculation(isUnrestricted);
    method_.setMolecularCharge(molecularCharge);
    method_.setSpinMultiplicity(spinMultiplicity);
    method_.setConvergenceCriteria(selfConsistenceCriterion);
    method_.setMaxIterations(maxIterations);
    method_.setScfMixer(scfMixerType);

    // After call to verifyPesValidity, update settings if they were internally changed due to
    // charge/spin multiplicity mismatch if not empty
    if (getLcaoMethod().getNumberAtomicOrbitals() != 0) {
      method_.verifyPesValidity();
      settings_->modifyBool(Utils::SettingsNames::unrestrictedCalculation, method_.unrestrictedCalculationRunning());
      settings_->modifyInt(Utils::SettingsNames::molecularCharge, method_.getMolecularCharge());
      settings_->modifyInt(Utils::SettingsNames::spinMultiplicity, method_.spinMultiplicity());
    }
  }
  else {
    throw Core::InitializationException("settings invalid!");
  }
}

std::string DFTB2MethodWrapper::name() const {
  return "DFTB2";
}

void DFTB2MethodWrapper::initialize() {
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

const Utils::LcaoMethod& DFTB2MethodWrapper::getLcaoMethod() const {
  return method_;
}

Utils::LcaoMethod& DFTB2MethodWrapper::getLcaoMethod() {
  return method_;
}

void DFTB2MethodWrapper::calculateImpl(Utils::derivativeType requiredDerivative) {
  method_.calculate(requiredDerivative);
}

Utils::DensityMatrix DFTB2MethodWrapper::getDensityMatrixGuess() const {
  return method_.getDensityMatrixGuess();
}

bool DFTB2MethodWrapper::successfulCalculation() const {
  return method_.hasConverged();
}

} /* namespace Sparrow */
} /* namespace Scine */
