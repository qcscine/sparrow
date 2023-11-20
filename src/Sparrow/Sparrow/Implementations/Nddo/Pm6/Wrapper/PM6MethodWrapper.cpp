/**
 * @file PM6MethodWrapper.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Internal Includes */
#include "PM6MethodWrapper.h"
#include "PM6Settings.h"
#include <Sparrow/Implementations/Nddo/Parameters.h>
#include <Sparrow/Implementations/Nddo/TimeDependent/LinearResponse/CISData.h>
#include <Sparrow/Implementations/Nddo/Utils/DipoleUtils/NDDODipoleMatrixCalculator.h>
#include <Sparrow/Implementations/Nddo/Utils/DipoleUtils/NDDODipoleMomentCalculator.h>
#include <Sparrow/Implementations/Nddo/Utils/NDDOInitializer.h>
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

PM6MethodWrapper::PM6MethodWrapper() {
  dipoleMatrixCalculator_ = NDDODipoleMatrixCalculator<nddo::PM6Method>::create(method_);
  dipoleCalculator_ = NDDODipoleMomentCalculator<nddo::PM6Method>::create(method_, *dipoleMatrixCalculator_);
  this->settings_ = std::make_unique<PM6Settings>();
  requiredProperties_ = Utils::Property::Energy;
  applySettings();
}

PM6MethodWrapper::PM6MethodWrapper(const PM6MethodWrapper& rhs) : PM6MethodWrapper() {
  copyInto(*this, rhs);
}

PM6MethodWrapper& PM6MethodWrapper::operator=(const PM6MethodWrapper& rhs) {
  copyInto(*this, rhs);
  return *this;
}

PM6MethodWrapper::~PM6MethodWrapper() = default;

void PM6MethodWrapper::applySettings() {
  // Set whether the NDDO dipole approximation is used.
  auto useNDDOApprox = settings_->getBool(Utils::SettingsNames::NDDODipoleApproximation);
  auto& NDDODipoleCalculator = dynamic_cast<NDDODipoleMomentCalculator<nddo::PM6Method>&>(*dipoleCalculator_);
  NDDODipoleCalculator.useNDDOApproximation(useNDDOApprox);

  NDDOMethodWrapper::applySettings(settings_, method_);
}

std::string PM6MethodWrapper::name() const {
  return "PM6";
}

void PM6MethodWrapper::initialize() {
  try {
    auto parameterFile = settings_->getString(Utils::SettingsNames::methodParameters);

    if (parameterFile.empty()) {
      method_.getInitializer().getRawParameters() = nddo::pm6();
    }
    else {
      method_.readParameters(parameterFile);
    }
    method_.initialize();
  }
  catch (Utils::Methods::InitializationException& e) {
    throw Core::InitializationException(e.what());
  }
}

const Utils::LcaoMethod& PM6MethodWrapper::getLcaoMethod() const {
  return method_;
}

Utils::LcaoMethod& PM6MethodWrapper::getLcaoMethod() {
  return method_;
}

void PM6MethodWrapper::calculateImpl(Utils::Derivative requiredDerivative) {
  method_.calculate(requiredDerivative, getLog());
}

Eigen::MatrixXd PM6MethodWrapper::getOneElectronMatrix() const {
  return method_.getOneElectronMatrix().getMatrix();
}

Utils::SpinAdaptedMatrix PM6MethodWrapper::getTwoElectronMatrix() const {
  Utils::SpinAdaptedMatrix twoElectronMatrix;
  if (!method_.unrestrictedCalculationRunning()) {
    twoElectronMatrix = Utils::SpinAdaptedMatrix::createRestricted(method_.getTwoElectronMatrix().getMatrix());
  }
  else {
    twoElectronMatrix = Utils::SpinAdaptedMatrix::createUnrestricted(method_.getTwoElectronMatrix().getAlpha(),
                                                                     method_.getTwoElectronMatrix().getBeta());
  }
  return twoElectronMatrix;
}

Utils::DensityMatrix PM6MethodWrapper::getDensityMatrixGuess() const {
  return method_.getDensityMatrixGuess();
}

CISData PM6MethodWrapper::getCISDataImpl() const {
  return CISData::constructCISDataFromNDDOMethod(method_);
}

void PM6MethodWrapper::addElectronicContribution(std::shared_ptr<Utils::AdditiveElectronicContribution> contribution) {
  method_.addElectronicContribution(std::move(contribution));
}

bool PM6MethodWrapper::successfulCalculation() const {
  return method_.hasConverged();
}

} /* namespace Sparrow */
} /* namespace Scine */
