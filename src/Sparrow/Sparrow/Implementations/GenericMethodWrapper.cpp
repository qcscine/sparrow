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
#include "MoldenFileGenerator.h"
#include <Utils/CalculatorBasics/PropertyList.h>
/* External Includes */
#include <Core/Exceptions.h>
#include <Sparrow/StatesHandling/SparrowState.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/GeometricDerivatives/NumericalHessianCalculator.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/Scf/MethodExceptions.h>
#include <Utils/Scf/MethodInterfaces/LcaoMethod.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <boost/filesystem.hpp>
#include <fstream>
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
  bool requiredHessian = requiredProperties_.containsSubSet(Utils::Property::Hessian) ||
                         requiredProperties_.containsSubSet(Utils::Property::Thermochemistry);
  bool requiredDipoleGradient = requiredProperties_.containsSubSet(Utils::Property::DipoleGradient);
  bool requiredSemiNumericalHessian = requiredHessian && !canCalculateAnalyticalHessian();

  results_ = Utils::Results{};
  if (dipoleMatrixCalculator_)
    dipoleMatrixCalculator_->invalidate();
  applySettings();
  calculateImpl(requiredDerivative);

  // If you want the Hessian, but cannot calculate it analytically,
  // calculate it semi-numerically. Same with Dipole Gradient.
  if (requiredSemiNumericalHessian || requiredDipoleGradient) {
    Utils::NumericalHessianCalculator hessianCalculator(*this);
    if (requiredDipoleGradient) {
      hessianCalculator.requiredDipoleGradient(true);
    }
    auto numericalResult = hessianCalculator.calculate();
    results_.set<Utils::Property::Hessian>(numericalResult.take<Utils::Property::Hessian>());
    if (requiredDipoleGradient) {
      results_.set<Utils::Property::DipoleGradient>(numericalResult.take<Utils::Property::DipoleGradient>());
    }
  }
  // If you want it, calculate and save the dipole.
  if (requiredProperties_.containsSubSet(Utils::Property::Dipole)) {
    results_.set<Utils::Property::Dipole>(dipoleCalculator_->calculate());
  }

  results_.set<Utils::Property::SuccessfulCalculation>(successfulCalculation());

  assembleResults(description);

  return results_;
};

void GenericMethodWrapper::setStructure(const Utils::AtomCollection& structure) {
  getLcaoMethod().setAtomCollection(structure);
  // applySettings() is here to prevent logging wrong settings due to the call to initialize.
  applySettings();
  initialize();
}

std::unique_ptr<Utils::AtomCollection> GenericMethodWrapper::getStructure() const {
  return std::make_unique<Utils::AtomCollection>(getLcaoMethod().getElementTypes(), getLcaoMethod().getPositions());
}

void GenericMethodWrapper::modifyPositions(Utils::PositionCollection newPositions) {
  if (getLcaoMethod().getElementTypes().size() != newPositions.rows()) {
    throw std::runtime_error("Position/ElementTypeCollection dimensionality mismatch.");
  }
  getLcaoMethod().setPositions(std::move(newPositions));
}

const Utils::PositionCollection& GenericMethodWrapper::getPositions() const {
  return getLcaoMethod().getPositions();
}

void GenericMethodWrapper::setRequiredProperties(const Utils::PropertyList& requiredProperties) {
  requiredProperties_ = requiredProperties;
}

Utils::PropertyList GenericMethodWrapper::getRequiredProperties() const {
  return requiredProperties_;
}

Utils::PropertyList GenericMethodWrapper::possibleProperties() const {
  return Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian |
         Utils::Property::BondOrderMatrix | Utils::Property::DensityMatrix | Utils::Property::AtomicCharges |
         Utils::Property::OverlapMatrix | Utils::Property::OrbitalEnergies | Utils::Property::CoefficientMatrix |
         Utils::Property::ElectronicOccupation | Utils::Property::Thermochemistry | Utils::Property::Description |
         Utils::Property::Dipole | Utils::Property::SuccessfulCalculation;
}

Utils::derivativeType GenericMethodWrapper::highestDerivativeRequired() const {
  // If gradient or hessian is contained in the requiredProperties_, then calculate the corresponding derivative.
  bool gradientsRequired = requiredProperties_.containsSubSet(Utils::Property::Gradients);
  bool hessianRequired = requiredProperties_.containsSubSet(Utils::Property::Hessian) ||
                         requiredProperties_.containsSubSet(Utils::Property::Thermochemistry);
  Utils::derivativeType requiredDerivative = Utils::derivativeType::none;
  if (gradientsRequired && possibleProperties().containsSubSet(Utils::Property::Gradients)) {
    requiredDerivative = Utils::derivativeType::first;
  }
  if (hessianRequired && canCalculateAnalyticalHessian()) {
    requiredDerivative = Utils::derivativeType::second_full;
  }
  return requiredDerivative;
}

void GenericMethodWrapper::assembleResults(const std::string& description) {
  results_.set<Utils::Property::Description>(description);
  results_.set<Utils::Property::Energy>(getLcaoMethod().getEnergy());

  if (requiredProperties_.containsSubSet(Utils::Property::Gradients)) {
    results_.set<Utils::Property::Gradients>(getLcaoMethod().getGradients());
  }

  if (highestDerivativeRequired() >= Utils::derivativeType::second_full) {
    results_.set<Utils::Property::Hessian>(getLcaoMethod().getFullSecondDerivatives().getHessianMatrix());
  }

  if (requiredProperties_.containsSubSet(Utils::Property::BondOrderMatrix)) {
    results_.set<Utils::Property::BondOrderMatrix>(getLcaoMethod().getBondOrderCollection());
  }

  if (requiredProperties_.containsSubSet(Utils::Property::DensityMatrix)) {
    results_.set<Utils::Property::DensityMatrix>(getLcaoMethod().getDensityMatrix());
  }

  if (requiredProperties_.containsSubSet(Utils::Property::AtomicCharges)) {
    results_.set<Utils::Property::AtomicCharges>(getLcaoMethod().getAtomicCharges());
  }

  if (requiredProperties_.containsSubSet(Utils::Property::OverlapMatrix)) {
    results_.set<Utils::Property::OverlapMatrix>(getLcaoMethod().getOverlapMatrix());
  }

  if (requiredProperties_.containsSubSet(Utils::Property::OrbitalEnergies)) {
    results_.set<Utils::Property::OrbitalEnergies>(getLcaoMethod().getSingleParticleEnergies());
  }

  if (requiredProperties_.containsSubSet(Utils::Property::CoefficientMatrix) ||
      requiredProperties_.containsSubSet(Utils::Property::Thermochemistry)) {
    results_.set<Utils::Property::CoefficientMatrix>(getLcaoMethod().getMolecularOrbitals());
  }

  if (requiredProperties_.containsSubSet(Utils::Property::ElectronicOccupation) ||
      requiredProperties_.containsSubSet(Utils::Property::Thermochemistry)) {
    results_.set<Utils::Property::ElectronicOccupation>(getLcaoMethod().getElectronicOccupation());
  }

  Utils::ResultsAutoCompleter resultsAutoCompleter(*getStructure());
  resultsAutoCompleter.setCoreCharges(getLcaoMethod().getAtomicCharges());
  auto zpveInclusion = getZPVEInclusion() ? Utils::ZPVEInclusion::alreadyIncluded : Utils::ZPVEInclusion::notIncluded;
  resultsAutoCompleter.setZPVEInclusion(zpveInclusion);

  if (requiredProperties_.containsSubSet(Utils::Property::Thermochemistry)) {
    resultsAutoCompleter.addOneWantedProperty(Utils::Property::Thermochemistry);
    resultsAutoCompleter.setTemperature(settings().getDouble(Utils::SettingsNames::temperature));
  }

  resultsAutoCompleter.generateProperties(results_, *getStructure());
}

bool GenericMethodWrapper::supportsMethodFamily(const std::string& methodFamily) const {
  return methodFamily == name();
}

void GenericMethodWrapper::loadState(std::shared_ptr<Core::State> state) {
  auto sparrowState = std::dynamic_pointer_cast<SparrowState>(state);
  if (!sparrowState)
    throw Core::StateCastingException();
  if (getLcaoMethod().getNumberElectrons() != sparrowState->getDensityMatrix().numberElectrons()) {
    throw IncompatibleStateException();
  }
  getLcaoMethod().setDensityMatrix(sparrowState->getDensityMatrix());
}

std::shared_ptr<Core::State> GenericMethodWrapper::getState() const {
  return std::make_shared<SparrowState>(getLcaoMethod().getDensityMatrix());
}

std::string GenericMethodWrapper::getStoNGExpansionPath() const {
  auto path = Utils::NativeFilenames::combinePathSegments(settings().getString(Utils::SettingsNames::parameterRootDirectory),
                                                          settings().getString(Utils::SettingsNames::parameterFile));
  boost::filesystem::path methodRoot(path);
  auto pathToBasis = methodRoot.parent_path();

  if (name() == "DFTB0" || name() == "DFTB2" || name() == "DFTB3")
    pathToBasis = pathToBasis.parent_path();

  pathToBasis = pathToBasis / "STO-6G.basis";
  return pathToBasis.string();
}

bool GenericMethodWrapper::canCalculateAnalyticalHessian() const {
  return false;
}

void GenericMethodWrapper::generateWavefunctionInformation(const std::string& filename) {
  std::ofstream fileOut(filename);
  if (!fileOut.is_open())
    throw std::runtime_error("Impossible to open file " + filename + ".");
  generateWavefunctionInformation(fileOut);
}

void GenericMethodWrapper::generateWavefunctionInformation(std::ostream& out) {
  getLcaoMethod().calculate(Utils::derivativeType::none);
  MoldenFileGenerator sparrow2molden(*this);
  sparrow2molden.generateWavefunctionInformation(out);
}

void GenericMethodWrapper::applySettings() {
}

} /* namespace Sparrow */
} // namespace Scine
