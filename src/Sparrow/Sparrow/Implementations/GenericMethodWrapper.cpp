/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Internal Includes */
#include "GenericMethodWrapper.h"
#include "DipoleMatrixCalculator.h"
#include "DipoleMomentCalculator.h"
/* External Includes */
#include <Core/Exceptions.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/GeometricDerivatives/NumericalHessianCalculator.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/MethodEssentials/Methods/LCAOMethod.h>
#include <Utils/MethodEssentials/util/MethodExceptions.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <memory>

namespace Scine {
namespace Sparrow {

GenericMethodWrapper::GenericMethodWrapper() = default;

GenericMethodWrapper::~GenericMethodWrapper() = default;

Utils::Settings& GenericMethodWrapper::settings() {
  return *settings_;
}

const Utils::Settings& GenericMethodWrapper::settings() const {
  return *settings_;
}

Utils::StatesHandler& GenericMethodWrapper::statesHandler() {
  return *statesHandler_;
}

const Utils::StatesHandler& GenericMethodWrapper::statesHandler() const {
  return *statesHandler_;
}

Utils::Results& GenericMethodWrapper::results() {
  return results_;
}

const Utils::Results& GenericMethodWrapper::results() const {
  return results_;
}

const Utils::Results& GenericMethodWrapper::calculate(std::string description) {
  if (getPositions().rows() == 0) {
    throw Core::EmptyMolecularStructureException();
  }
  auto requiredDerivative = highestDerivativeRequired();
  bool requiredDipoleGradient = requiredProperties_.containsSubSet(Utils::Property::DipoleGradient);
  bool requiredSemiNumericalHessian = requiredProperties_.containsSubSet(Utils::Property::Hessian) &&
                                      !this->possibleProperties().containsSubSet(Utils::Property::Hessian);

  if (dipoleMatrixCalculator_)
    dipoleMatrixCalculator_->invalidate();
  applySettings();
  calculateImpl(requiredDerivative);

  results_ = assembleResults(description);

  // If you want the Hessian, but cannot calculate it analytically,
  // calculate it semi-numerically. Same with Dipole Gradient.
  if (requiredSemiNumericalHessian || requiredDipoleGradient) {
    Utils::NumericalHessianCalculator hessianCalculator(*this);
    if (requiredDipoleGradient) {
      hessianCalculator.requiredDipoleGradient(true);
    }
    auto numericalResult = hessianCalculator.calculate();
    results_.setHessian(numericalResult.takeHessian());
    if (requiredDipoleGradient) {
      results_.setDipoleGradient(numericalResult.takeDipoleGradient());
    }
  }
  // If you want it, calculate and save the dipole.
  if (requiredProperties_.containsSubSet(Utils::Property::Dipole)) {
    results_.setDipole(dipoleCalculator_->calculate());
  }

  return results_;
};

void GenericMethodWrapper::setStructure(const Utils::AtomCollection& structure) {
  getLCAOMethod().setAtomCollection(structure);
  // applySettings() is here to prevent logging wrong settings due to the call to initialize.
  applySettings();
  initialize();
}

std::unique_ptr<Utils::AtomCollection> GenericMethodWrapper::getStructure() const {
  return std::make_unique<Utils::AtomCollection>(getLCAOMethod().getElementTypes(), getLCAOMethod().getPositions());
}

void GenericMethodWrapper::modifyPositions(Utils::PositionCollection newPositions) {
  if (getLCAOMethod().getElementTypes().size() != newPositions.rows()) {
    throw std::runtime_error("Position/ElementTypeCollection dimensionality mismatch.");
  }
  getLCAOMethod().setPositions(std::move(newPositions));
}

const Utils::PositionCollection& GenericMethodWrapper::getPositions() const {
  return getLCAOMethod().getPositions();
}

void GenericMethodWrapper::setRequiredProperties(const Utils::PropertyList& requiredProperties) {
  requiredProperties_ = requiredProperties;
}

Utils::PropertyList GenericMethodWrapper::possibleProperties() const {
  return Utils::Property::Energy | Utils::Property::Gradients;
}

Utils::derivativeType GenericMethodWrapper::highestDerivativeRequired() const {
  // If gradient or hessian is contained in the requiredProperties_, then calculate the corresponding derivative.
  bool gradientsRequired = requiredProperties_.containsSubSet(Utils::Property::Gradients);
  bool hessianRequired = requiredProperties_.containsSubSet(Utils::Property::Hessian);
  Utils::derivativeType requiredDerivative = Utils::derivativeType::none;
  if (gradientsRequired && possibleProperties().containsSubSet(Utils::Property::Gradients)) {
    requiredDerivative = Utils::derivativeType::first;
  }
  if (hessianRequired && possibleProperties().containsSubSet(Utils::Property::Hessian)) {
    requiredDerivative = Utils::derivativeType::second_full;
  }
  return requiredDerivative;
}

Utils::Results GenericMethodWrapper::assembleResults(const std::string& description) const {
  Utils::Results results;
  results.setDescription(description);
  results.setEnergy(getLCAOMethod().getEnergy());

  if (requiredProperties_.containsSubSet(Utils::Property::Gradients)) {
    results.setGradients(getLCAOMethod().getGradients());
  }

  if (highestDerivativeRequired() >= Utils::derivativeType::second_full) {
    results.setHessian(getLCAOMethod().getFullSecondDerivatives().getHessianMatrix());
  }

  if (requiredProperties_.containsSubSet(Utils::Property::BondOrderMatrix)) {
    results.setBondOrders(getLCAOMethod().getBondOrderCollection());
  }

  return results;
}
} /* namespace Sparrow */
} // namespace Scine
