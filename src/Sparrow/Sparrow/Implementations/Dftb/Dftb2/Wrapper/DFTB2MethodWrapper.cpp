/**
 * @file DFTB2MethodWrapper.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Internal Includes */
#include "DFTB2MethodWrapper.h"
#include "DFTB2Settings.h"
#include <Sparrow/Implementations/Dftb/TimeDependent/LinearResponse/TDDFTBData.h>
#include <Sparrow/Implementations/Dftb/Utils/DipoleUtils/DFTBDipoleMatrixCalculator.h>
#include <Sparrow/Implementations/Dftb/Utils/DipoleUtils/DFTBDipoleMomentCalculator.h>
/* External Includes */
#include <Core/Exceptions.h>
#include <Sparrow/Implementations/Dftb/Utils/SecondOrderFock.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <Utils/Scf/MethodExceptions.h>
#include <Utils/Scf/MethodInterfaces/AdditiveElectronicContribution.h>
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
  if (settings_->valid()) {
    auto spinMode = Utils::SpinModeInterpreter::getSpinModeFromString(settings_->getString(Utils::SettingsNames::spinMode));
    int molecularCharge = settings_->getInt(Utils::SettingsNames::molecularCharge);
    int spinMultiplicity = settings_->getInt(Utils::SettingsNames::spinMultiplicity);
    double selfConsistenceCriterion = settings_->getDouble(Utils::SettingsNames::selfConsistenceCriterion);
    double densityRmsdThreshold = settings_->getDouble(Utils::SettingsNames::densityRmsdCriterion);
    int maxScfIterations = settings_->getInt(Utils::SettingsNames::maxScfIterations);
    auto scfMixerName = settings_->getString(Utils::SettingsNames::mixer);
    auto scfMixerType = Utils::UniversalSettings::SettingPopulator::stringToScfMixer(scfMixerName);

    if (spinMode == Utils::SpinMode::Any) {
      if (spinMultiplicity == 1) {
        method_.setUnrestrictedCalculation(false);
      }
      else {
        method_.setUnrestrictedCalculation(true);
      }
    }
    if (spinMode == Utils::SpinMode::Restricted) {
      if (spinMultiplicity != 1) {
        throw std::runtime_error("Restricted calculation with multiplicity != 1 requested.");
      }
      method_.setUnrestrictedCalculation(false);
    }
    if (spinMode == Utils::SpinMode::Unrestricted) {
      method_.setUnrestrictedCalculation(true);
    }
    method_.setMolecularCharge(molecularCharge);
    method_.setSpinMultiplicity(spinMultiplicity);
    method_.setConvergenceCriteria({selfConsistenceCriterion, densityRmsdThreshold});
    method_.setMaxIterations(maxScfIterations);
    method_.setScfMixer(scfMixerType);
  }
  else {
    settings_->throwIncorrectSettings();
  }
}

std::string DFTB2MethodWrapper::name() const {
  return "DFTB2";
}

void DFTB2MethodWrapper::initialize() {
  try {
    auto parameters = settings_->getString(Utils::SettingsNames::methodParameters);
    method_.initializeFromParameterPath(parameters);
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

void DFTB2MethodWrapper::calculateImpl(Utils::Derivative requiredDerivative) {
  method_.calculate(requiredDerivative, getLog());
}

Utils::DensityMatrix DFTB2MethodWrapper::getDensityMatrixGuess() const {
  return method_.getDensityMatrixGuess();
}

void DFTB2MethodWrapper::addElectronicContribution(std::shared_ptr<Utils::AdditiveElectronicContribution> contribution) {
  method_.addElectronicContribution(std::move(contribution));
}

TDDFTBData DFTB2MethodWrapper::getTDDFTBDataImpl() const {
  return TDDFTBData::constructTDDFTBDataFromDFTBMethod(method_);
}

bool DFTB2MethodWrapper::successfulCalculation() const {
  return method_.hasConverged();
}

} /* namespace Sparrow */
} /* namespace Scine */
