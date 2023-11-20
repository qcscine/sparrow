/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Internal Includes */
#include "GenericMethodWrapper.h"
#include "DipoleMatrixCalculator.h"
#include "DipoleMomentCalculator.h"
#include "MoldenFileGenerator.h"
#include "Sto6gParameters.h"
#include "Utils/Scf/LcaoUtils/SpinMode.h"
#include <Utils/CalculatorBasics/PropertyList.h>
/* External Includes */
#include <Core/Exceptions.h>
#include <Sparrow/StatesHandling/SparrowState.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/GeometricDerivatives/NumericalHessianCalculator.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/ElementInfo.h>
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

GenericMethodWrapper::GenericMethodWrapper(const GenericMethodWrapper& rhs) : CloneInterface(rhs) {
}

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
  // Check method and basis set fields
  checkBasicSettings();
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
}

void GenericMethodWrapper::setStructure(const Utils::AtomCollection& structure) {
  results_ = {};
  getLcaoMethod().setAtomCollection(structure);
  // Apply the settings used in intialization, as molecular charge
  applySettings();
  initialize();
}

std::unique_ptr<Utils::AtomCollection> GenericMethodWrapper::getStructure() const {
  return std::make_unique<Utils::AtomCollection>(getLcaoMethod().getElementTypes(), getLcaoMethod().getPositions());
}

void GenericMethodWrapper::modifyPositions(Utils::PositionCollection newPositions) {
  results_ = {};
  if (int(getLcaoMethod().getElementTypes().size()) != int(newPositions.rows())) {
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
         Utils::Property::AtomicHessians | Utils::Property::BondOrderMatrix | Utils::Property::DensityMatrix |
         Utils::Property::AtomicCharges | Utils::Property::OverlapMatrix | Utils::Property::OrbitalEnergies |
         Utils::Property::CoefficientMatrix | Utils::Property::ElectronicOccupation | Utils::Property::Thermochemistry |
         Utils::Property::Description | Utils::Property::Dipole | Utils::Property::SuccessfulCalculation;
}

Utils::Derivative GenericMethodWrapper::highestDerivativeRequired() const {
  // If gradient or hessian is contained in the requiredProperties_, then calculate the corresponding derivative.
  bool gradientsRequired = requiredProperties_.containsSubSet(Utils::Property::Gradients);
  bool hessianRequired = requiredProperties_.containsSubSet(Utils::Property::Hessian) ||
                         requiredProperties_.containsSubSet(Utils::Property::Thermochemistry);
  bool atomicHessiansRequired = requiredProperties_.containsSubSet(Utils::Property::AtomicHessians);
  Utils::Derivative requiredDerivative = Utils::Derivative::None;
  if (gradientsRequired && possibleProperties().containsSubSet(Utils::Property::Gradients)) {
    requiredDerivative = Utils::Derivative::First;
  }
  if (hessianRequired && canCalculateAnalyticalHessian()) {
    requiredDerivative = Utils::Derivative::SecondFull;
  }
  if (atomicHessiansRequired) {
    requiredDerivative = Utils::Derivative::SecondAtomic;
  }
  return requiredDerivative;
}

void GenericMethodWrapper::checkBasicSettings() {
  Utils::Settings settingsCopy(*(this->settings_));
  settingsCopy.resetToDefaults();
  auto methodDefault = settingsCopy.getString(Utils::SettingsNames::method);
  auto currentMethod = this->settings_->getString(Utils::SettingsNames::method);

  // Check selected method
  std::transform(currentMethod.begin(), currentMethod.end(), currentMethod.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  if (currentMethod != methodDefault && currentMethod != "any") {
    throw std::runtime_error("This calculator does not provide the requested method.");
  }
  settings_->modifyString(Utils::SettingsNames::method, methodDefault);
}

void GenericMethodWrapper::assembleResults(const std::string& description) {
  results_.set<Utils::Property::Description>(description);
  results_.set<Utils::Property::Energy>(getLcaoMethod().getEnergy());

  if (requiredProperties_.containsSubSet(Utils::Property::Gradients)) {
    Utils::GradientCollection grad(getPositions().rows(), 3);
    if (highestDerivativeRequired() == Utils::Derivative::First) {
      grad = getLcaoMethod().getGradients();
    }
    else if (highestDerivativeRequired() == Utils::Derivative::SecondAtomic) {
      int index = 0;
      for (const auto& sd : getLcaoMethod().getAtomicSecondDerivatives()) {
        grad.row(index++) = sd.deriv().transpose();
      }
    }
    else if (highestDerivativeRequired() == Utils::Derivative::SecondFull) {
      grad = getLcaoMethod().getFullSecondDerivatives().getReferenceGradients();
    }
    results_.set<Utils::Property::Gradients>(std::move(grad));
  }

  if (highestDerivativeRequired() >= Utils::Derivative::SecondFull) {
    results_.set<Utils::Property::Hessian>(getLcaoMethod().getFullSecondDerivatives().getHessianMatrix());
  }

  if (requiredProperties_.containsSubSet(Utils::Property::AtomicHessians)) {
    results_.set<Utils::Property::AtomicHessians>(getLcaoMethod().getAtomicSecondDerivatives());
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
    results_.set<Utils::Property::OverlapMatrix>(getLcaoMethod().getOverlapMatrix().selfadjointView<Eigen::Lower>());
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

  if (requiredProperties_.containsSubSet(Utils::Property::AOtoAtomMapping)) {
    results_.set<Utils::Property::AOtoAtomMapping>(getLcaoMethod().getAtomsOrbitalsIndexesHolder());
  }

  if (requiredProperties_.containsSubSet(Utils::Property::AtomicGtos)) {
    results_.set<Utils::Property::AtomicGtos>(getAtomicGtosMap());
  }

  Utils::ResultsAutoCompleter resultsAutoCompleter(*getStructure());
  resultsAutoCompleter.setCoreCharges(getLcaoMethod().getAtomicCharges());
  auto zpveInclusion = getZPVEInclusion() ? Utils::ZPVEInclusion::alreadyIncluded : Utils::ZPVEInclusion::notIncluded;
  resultsAutoCompleter.setZPVEInclusion(zpveInclusion);

  if (requiredProperties_.containsSubSet(Utils::Property::Thermochemistry)) {
    resultsAutoCompleter.addOneWantedProperty(Utils::Property::Thermochemistry);
    resultsAutoCompleter.setTemperature(settings().getDouble(Utils::SettingsNames::temperature));
    resultsAutoCompleter.setPressure(settings().getDouble(Utils::SettingsNames::pressure));
    resultsAutoCompleter.setMolecularSymmetryNumber(settings().getInt(Utils::SettingsNames::symmetryNumber));
  }

  resultsAutoCompleter.generateProperties(results_, *getStructure());

  // update spin_mode
  if (settings_->valueExists(Utils::SettingsNames::spinMode)) {
    Utils::SpinMode spinmode =
        getLcaoMethod().unrestrictedCalculationRunning() ? Utils::SpinMode::Unrestricted : Utils::SpinMode::Restricted;
    settings_->modifyString(Utils::SettingsNames::spinMode, Utils::SpinModeInterpreter::getStringFromSpinMode(spinmode));
  }
}

bool GenericMethodWrapper::supportsMethodFamily(const std::string& methodFamily) const {
  return methodFamily == name();
}

void GenericMethodWrapper::loadState(std::shared_ptr<Core::State> state) {
  auto sparrowState = std::dynamic_pointer_cast<SparrowState>(state);
  if (!sparrowState)
    throw Core::StateCastingException();
  if (sparrowState->getDensityMatrix().numberElectrons() != 0) {
    if (getLcaoMethod().getNumberElectrons() != sparrowState->getDensityMatrix().numberElectrons()) {
      throw IncompatibleStateException();
    }
    getLcaoMethod().setDensityMatrix(sparrowState->getDensityMatrix());
  }
}

std::shared_ptr<Core::State> GenericMethodWrapper::getState() const {
  return std::make_shared<SparrowState>(getLcaoMethod().getDensityMatrix());
}

std::string GenericMethodWrapper::getStoNGExpansionPath() const {
  boost::filesystem::path methodRoot(settings().getString(Utils::SettingsNames::methodParameters));
  auto pathToBasis = methodRoot.parent_path();

  if (name() == "DFTB0" || name() == "DFTB2" || name() == "DFTB3") {
    pathToBasis = pathToBasis.parent_path();
    pathToBasis = pathToBasis / "STO-6G.basis";
  }
  else {
    std::string fileName = name() + "-STO-6G.basis";
    pathToBasis = pathToBasis / fileName;
  }
  return pathToBasis.string();
}

std::unordered_map<int, Utils::AtomicGtos> GenericMethodWrapper::getAtomicGtosMap() const {
  std::unordered_map<int, Utils::AtomicGtos> expansionMap;
  if (name() == "MNDO") {
    expansionMap = Sto6g::mndo();
  }
  else if (name() == "AM1") {
    expansionMap = Sto6g::am1();
  }
  else if (name() == "RM1") {
    expansionMap = Sto6g::rm1();
  }
  else if (name() == "PM3") {
    expansionMap = Sto6g::pm3();
  }
  else if (name() == "PM6") {
    expansionMap = Sto6g::pm6();
  }
  else if (name().find("DFTB") != std::string::npos) {
    expansionMap = Sto6g::dftb();
  }
  else {
    throw std::runtime_error("Atomic GTO Property is not possible for method " + name() + ".");
  }
  /* expansionMap includes STO for all elements, now give out only for exisiting elements in calc */
  std::unordered_map<int, Utils::AtomicGtos> necessaryExpansionMap;
  Utils::ElementTypeCollection elements = getLcaoMethod().getElementTypes();
  /* remove duplicate elements */
  auto last = std::unique(elements.begin(), elements.end());
  elements.erase(last, elements.end());
  /* find elements and save STOs in map */
  for (const auto& ele : elements) {
    auto it = expansionMap.find(Utils::ElementInfo::Z(ele));
    if (it == expansionMap.end()) {
      throw std::runtime_error("Could not find GTO for element " + Utils::ElementInfo::symbol(ele) + ".");
    }
    necessaryExpansionMap.insert(*it);
  }
  return necessaryExpansionMap;
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
  if (!results().has<Utils::Property::SuccessfulCalculation>() || !results().get<Utils::Property::SuccessfulCalculation>()) {
    getLcaoMethod().calculate(Utils::Derivative::None, getLog());
  }
  MoldenFileGenerator sparrow2molden(*this);
  sparrow2molden.generateWavefunctionInformation(out);
}

void GenericMethodWrapper::applySettings() {
}

} /* namespace Sparrow */
} // namespace Scine
