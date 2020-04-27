/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Internal Includes */
#include "AM1TypeMethodWrapper.h"
#include "AM1Settings.h"
#include <Sparrow/Implementations/Nddo/Utils/DipoleUtils/NDDODipoleMatrixCalculator.h>
#include <Sparrow/Implementations/Nddo/Utils/DipoleUtils/NDDODipoleMomentCalculator.h>
#include <Sparrow/Implementations/Nddo/Utils/OneElectronMatrix.h>
#include <Sparrow/Implementations/Nddo/Utils/TwoElectronMatrix.h>
/* External Includes */
#include <Core/Exceptions.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/Scf/MethodExceptions.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <memory>

namespace Scine {
namespace Sparrow {

template<class AM1Type>
AM1TypeMethodWrapper<AM1Type>::AM1TypeMethodWrapper() {
  this->dipoleMatrixCalculator_ = NDDODipoleMatrixCalculator<nddo::AM1Method>::create(method_);
  this->dipoleCalculator_ = NDDODipoleMomentCalculator<nddo::AM1Method>::create(method_, *(this->dipoleMatrixCalculator_));
}

template<class AM1Type>
void AM1TypeMethodWrapper<AM1Type>::applySettings() {
  // Set whether the NDDO dipole approximation is used.
  auto useNDDOApprox = this->settings_->getBool(Utils::SettingsNames::NDDODipoleApproximation);
  auto& NDDODipoleCalculator = dynamic_cast<NDDODipoleMomentCalculator<nddo::AM1Method>&>(*(this->dipoleCalculator_));
  NDDODipoleCalculator.useNDDOApproximation(useNDDOApprox);

  auto& derived = static_cast<AM1Type&>(*this);
  NDDOMethodWrapper::applySettings(derived.settings_, derived.method_);
};

template<class AM1Type>
std::string AM1TypeMethodWrapper<AM1Type>::name() const {
  return AM1Type::model;
}

template<class AM1Type>
void AM1TypeMethodWrapper<AM1Type>::initialize() {
  auto& derived = static_cast<AM1Type&>(*this);
  try {
    auto parameterFile = derived.settings_->getString(Utils::SettingsNames::parameterFile);
    auto resourceDirectory = derived.settings_->getString(Utils::SettingsNames::parameterRootDirectory);
    auto fullPathToParameters = Utils::NativeFilenames::combinePathSegments(resourceDirectory, parameterFile);
    method_.readParameters(fullPathToParameters);
    method_.initialize();
  }
  catch (Utils::Methods::InitializationException& e) {
    throw Core::InitializationException(e.what());
  }
}

template<class AM1Type>
Utils::LcaoMethod& AM1TypeMethodWrapper<AM1Type>::getLcaoMethod() {
  return method_;
}

template<class AM1Type>
const Utils::LcaoMethod& AM1TypeMethodWrapper<AM1Type>::getLcaoMethod() const {
  return method_;
}

template<class AM1Type>
void AM1TypeMethodWrapper<AM1Type>::calculateImpl(Utils::derivativeType requiredDerivative) {
  method_.calculate(requiredDerivative);
}

template<class AM1Type>
Eigen::MatrixXd AM1TypeMethodWrapper<AM1Type>::getOneElectronMatrix() const {
  return method_.getOneElectronMatrix().getMatrix();
}

template<class AM1Type>
Utils::SpinAdaptedMatrix AM1TypeMethodWrapper<AM1Type>::getTwoElectronMatrix() const {
  Utils::SpinAdaptedMatrix twoElectronMatrix;
  if (method_.unrestrictedCalculationRunning()) {
    twoElectronMatrix = Utils::SpinAdaptedMatrix::createRestricted(method_.getTwoElectronMatrix().getMatrix());
  }
  else {
    twoElectronMatrix = Utils::SpinAdaptedMatrix::createUnrestricted(method_.getTwoElectronMatrix().getAlpha(),
                                                                     method_.getTwoElectronMatrix().getBeta());
  }
  return twoElectronMatrix;
}

template<class AM1Type>
Utils::DensityMatrix AM1TypeMethodWrapper<AM1Type>::getDensityMatrixGuess() const {
  return method_.getDensityMatrixGuess();
}

template<class AM1Type>
bool AM1TypeMethodWrapper<AM1Type>::successfulCalculation() const {
  return method_.hasConverged();
}

AM1MethodWrapper::AM1MethodWrapper() {
  requiredProperties_ = Utils::Property::Energy;
  this->settings_ = std::make_unique<AM1Settings>();
  applySettings();
}

AM1MethodWrapper::AM1MethodWrapper(const AM1MethodWrapper& rhs) : AM1MethodWrapper() {
  copyInto(*this, rhs);
}

AM1MethodWrapper& AM1MethodWrapper::operator=(const AM1MethodWrapper& rhs) {
  copyInto(*this, rhs);
  return *this;
}

RM1MethodWrapper::RM1MethodWrapper() {
  requiredProperties_ = Utils::Property::Energy;
  this->settings_ = std::make_unique<RM1Settings>();
  applySettings();
}

RM1MethodWrapper::RM1MethodWrapper(const RM1MethodWrapper& rhs) : RM1MethodWrapper() {
  copyInto(*this, rhs);
}

RM1MethodWrapper& RM1MethodWrapper::operator=(const RM1MethodWrapper& rhs) {
  copyInto(*this, rhs);
  return *this;
}

PM3MethodWrapper::PM3MethodWrapper() {
  requiredProperties_ = Utils::Property::Energy;
  this->settings_ = std::make_unique<PM3Settings>();
  applySettings();
}

PM3MethodWrapper::PM3MethodWrapper(const PM3MethodWrapper& rhs) : PM3MethodWrapper() {
  copyInto(*this, rhs);
}

PM3MethodWrapper& PM3MethodWrapper::operator=(const PM3MethodWrapper& rhs) {
  copyInto(*this, rhs);
  return *this;
}

AM1MethodWrapper::~AM1MethodWrapper() = default;
RM1MethodWrapper::~RM1MethodWrapper() = default;
PM3MethodWrapper::~PM3MethodWrapper() = default;

template<class AM1Type>
AM1TypeMethodWrapper<AM1Type>::~AM1TypeMethodWrapper() = default;

} /* namespace Sparrow */
} /* namespace Scine */
