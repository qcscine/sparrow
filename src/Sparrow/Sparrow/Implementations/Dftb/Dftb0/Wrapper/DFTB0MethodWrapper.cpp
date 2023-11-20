/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Internal Includes */
#include "DFTB0MethodWrapper.h"
#include "DFTB0Settings.h"
#include <Sparrow/Implementations/Dftb/TimeDependent/LinearResponse/TDDFTBData.h>
#include <Sparrow/Implementations/Dftb/Utils/DipoleUtils/DFTBDipoleMatrixCalculator.h>
#include <Sparrow/Implementations/Dftb/Utils/DipoleUtils/DFTBDipoleMomentCalculator.h>
#include <Sparrow/StatesHandling/SparrowState.h>
/* External Includes */
#include <Core/Exceptions.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/Scf/MethodExceptions.h>
#include <Utils/Scf/MethodInterfaces/AdditiveElectronicContribution.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <memory>

namespace Scine {
namespace Sparrow {

DFTB0MethodWrapper::DFTB0MethodWrapper() {
  this->settings_ = std::make_unique<DFTB0Settings>();
  requiredProperties_ = Utils::Property::Energy;
  dipoleCalculator_ = std::make_unique<DFTBDipoleMomentCalculator<dftb::DFTB0>>(method_);
  dipoleMatrixCalculator_ = DFTBDipoleMatrixCalculator<dftb::DFTB0>::create(method_);
  applySettings();
}

DFTB0MethodWrapper::DFTB0MethodWrapper(const DFTB0MethodWrapper& rhs) : DFTB0MethodWrapper() {
  copyInto(*this, rhs);
}

DFTB0MethodWrapper& DFTB0MethodWrapper::operator=(const DFTB0MethodWrapper& rhs) {
  copyInto(*this, rhs);
  return *this;
}

DFTB0MethodWrapper::~DFTB0MethodWrapper() = default;

void DFTB0MethodWrapper::applySettings() {
  GenericMethodWrapper::applySettings();
  if (settings_->valid()) {
    int molecularCharge = settings_->getInt(Utils::SettingsNames::molecularCharge);
    method_.setMolecularCharge(molecularCharge);
  }
  else {
    settings_->throwIncorrectSettings();
  }
}

std::string DFTB0MethodWrapper::name() const {
  return "DFTB0";
}

void DFTB0MethodWrapper::initialize() {
  try {
    const auto parameters = settings_->getString(Utils::SettingsNames::methodParameters);
    method_.initializeFromParameterPath(parameters);
  }
  catch (Utils::Methods::InitializationException& e) {
    throw Core::InitializationException(e.what());
  }
}

Utils::LcaoMethod& DFTB0MethodWrapper::getLcaoMethod() {
  return method_;
}

const Utils::LcaoMethod& DFTB0MethodWrapper::getLcaoMethod() const {
  return method_;
}

void DFTB0MethodWrapper::calculateImpl(Utils::Derivative requiredDerivative) {
  method_.calculate(requiredDerivative, getLog());
}

Utils::DensityMatrix DFTB0MethodWrapper::getDensityMatrixGuess() const {
  Utils::DensityMatrix densityMatrix;
  densityMatrix.resize(method_.getNumberAtomicOrbitals());
  return densityMatrix;
}

void DFTB0MethodWrapper::copyInto(DFTB0MethodWrapper& instance, const DFTB0MethodWrapper& classToCopy) {
  auto results = classToCopy.results();
  instance.settings() = classToCopy.settings();
  instance.setStructure(*classToCopy.getStructure());
  instance.results() = results;
}

void DFTB0MethodWrapper::addElectronicContribution(std::shared_ptr<Utils::AdditiveElectronicContribution> contribution) {
  method_.addElectronicContribution(std::move(contribution));
}

TDDFTBData DFTB0MethodWrapper::getTDDFTBDataImpl() const {
  return TDDFTBData::constructTDDFTBDataFromDFTBMethod(method_);
}

bool DFTB0MethodWrapper::successfulCalculation() const {
  return true;
}

void DFTB0MethodWrapper::loadState(std::shared_ptr<Core::State> /*state*/) {
}

} /* namespace Sparrow */
} /* namespace Scine */
