/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "UvVisCalculator.h"
#include "../SpectroscopySettings.h"
#include "../Utils/Spectrum.h"
#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/CalculatorWithReference.h>
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Sparrow/Implementations/RealTimeSpectroscopy/Utils/LineWidthGenerator.h>
#include <Sparrow/Implementations/TimeDependent/GuessPropagator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/TimeDependent/TransitionDipoleCalculator.h>
#include <Utils/UniversalSettings/SettingsNames.h>
namespace Scine {
namespace Sparrow {
namespace RealTimeSpectroscopy {

UvVisCalculator::UvVisCalculator() : settings_(std::make_shared<UvVisSettings>()) {
}

UvVisCalculator::UvVisCalculator(std::shared_ptr<Utils::Settings> settings) {
  settings_ = std::move(settings);
}

void UvVisCalculator::initialize(const Utils::ElementTypeCollection& elements) {
  createExcitedStatesCalculator();
  createReferenceCalculator();

  calculator_->getReferenceCalculator().setRequiredProperties(Utils::Property::Energy | Utils::Property::Dipole |
                                                              Utils::Property::DipoleMatrixAO);
  // Has to be accurate enough.
  calculator_->getReferenceCalculator().settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-9);
  if (!settings_->getString(Utils::SettingsNames::methodParameters).empty()) {
    calculator_->getReferenceCalculator().settings().modifyString(
        Utils::SettingsNames::methodParameters, settings_->getString(Utils::SettingsNames::methodParameters));
  }
  Utils::PositionCollection positions = Utils::PositionCollection::Zero(elements.size(), 3);
  auto structure = Utils::AtomCollection{elements, positions};
  calculator_->getReferenceCalculator().setStructure(structure);

  for (const auto& key : calculator_->settings().getKeys()) {
    if (settings_->valueExists(key)) {
      calculator_->settings().modifyValue(key, settings_->getValue(key));
    }
  }
  calculator_->applySettings();
  lastExcitedStates_.reset();
  guessPropagator_ = std::make_unique<GuessPropagator>();
  guessPropagator_->setDimension(settings_->getInt(uvVisGuessPropagatorDiisDimension));
}

void UvVisCalculator::createExcitedStatesCalculator() {
  auto m = settings_->getString(Utils::SettingsNames::method);
  std::string excitedStatesMethod;
  if (m == "PM6" || m == "MNDO" || m == "AM1" || m == "RM1" || m == "PM3")
    excitedStatesMethod = "CIS-NDDO";
  else if (m == "DFTB0" || m == "DFTB2" || m == "DFTB3")
    excitedStatesMethod = "TD-DFTB";
  else
    throw std::runtime_error("Method " + excitedStatesMethod + " not available.");

  auto& manager = Core::ModuleManager::getInstance();
  try {
    calculator_ =
        std::dynamic_pointer_cast<LinearResponseCalculator>(manager.get<Core::CalculatorWithReference>(excitedStatesMethod));
  }
  catch (std::exception& e) {
    throw std::runtime_error("Excited states method " + excitedStatesMethod + " not found.");
  }
}

void UvVisCalculator::createReferenceCalculator() {
  auto& manager = Core::ModuleManager::getInstance();
  try {
    calculator_->setReferenceCalculator(manager.get<Core::Calculator>(settings_->getString(Utils::SettingsNames::method)));
  }
  catch (std::exception& e) {
    throw std::runtime_error("Method " + settings_->getString(Utils::SettingsNames::method) + " not found.");
  }
  // calculator_->getLog().debug.add("cout", Core::Log::coutSink());
  calculator_->setLog(Core::Log::silent());
  calculator_->getLog().output.add("verbose", Core::Log::fileSink("verbose.out"));
  calculator_->getLog().warning.add("warn", Core::Log::cerrSink());
  calculator_->getLog().error.add("err", Core::Log::cerrSink());
}

void UvVisCalculator::modifyPositions(const Utils::PositionCollection& positions) {
  if (calculator_)
    calculator_->getReferenceCalculator().modifyPositions(positions);
}

void UvVisCalculator::updateState(std::shared_ptr<Core::State> state) {
  if (calculator_)
    calculator_->getReferenceCalculator().loadState(std::move(state));
}

Spectrum UvVisCalculator::calculate(const Utils::PositionCollection& positions) {
  calculator_->applySettings();

  if (settings_->getString("prune_basis") == "energy" && !sizeOutput_) {
    sizeOutput_ = std::make_unique<std::ofstream>("uvvis_dimensions.out");
  }

  calculator_->getReferenceCalculator().modifyPositions(positions);
  if (!settings_->getString(Utils::SettingsNames::methodParameters).empty()) {
    calculator_->getReferenceCalculator().settings().modifyString(
        Utils::SettingsNames::methodParameters, settings_->getString(Utils::SettingsNames::methodParameters));
  }
  calculator_->referenceCalculation();
  if (lastExcitedStates_) {
    auto guess = guessPropagator_->calculateGuessAtNewPosition(calculator_->getReferenceCalculator().getPositions());
    calculator_->setGuess(guess);
  }
  auto results = calculator_->calculate();
  std::pair<std::string, std::string> labels = {"Vertical Transition Energy / eV",
                                                "Integral Absorption Coefficient / km mol^-1"};

  if (!results.has<Utils::Property::ExcitedStates>()) {
    throw std::runtime_error("No excited states in the results!");
  }
  const auto& excitedStatesResults = results.get<Utils::Property::ExcitedStates>();
  std::shared_ptr<Utils::ElectronicTransitionResult> spinResolvedResults;
  if (excitedStatesResults.singlet)
    spinResolvedResults = excitedStatesResults.singlet;
  else if (excitedStatesResults.unrestricted)
    spinResolvedResults = excitedStatesResults.unrestricted;
  else
    throw std::runtime_error("No singlet or unrestricted excited states in the results!");

  if (settings_->getString("prune_basis") == "energy" && sizeOutput_) {
    *sizeOutput_ << spinResolvedResults->eigenStates.eigenVectors.rows() << std::endl;
  }

  Eigen::VectorXd energies = spinResolvedResults->eigenStates.eigenValues;
  Eigen::VectorXd intensities = Utils::TransitionDipoleCalculator::transitionDipoleMomentToOscillatorStrength(
      spinResolvedResults->transitionDipoles, energies);
  lastExcitedStates_ = calculator_->getGuess();
  guessPropagator_->record(calculator_->getReferenceCalculator().getPositions());
  guessPropagator_->record(*lastExcitedStates_);
  Spectrum spectrum{energies * Utils::Constants::ev_per_hartree, intensities, labels};
  LineWidthGenerator processor(spectrum);
  return processor.generateLorentzianProfile(settings_->getDouble(resolutionOption), settings_->getDouble(fwhmOption));
}

Utils::Settings& UvVisCalculator::settings() {
  return *settings_;
}

const Utils::Settings& UvVisCalculator::settings() const {
  return *settings_;
}

UvVisCalculator::~UvVisCalculator() = default;

} // namespace RealTimeSpectroscopy
} // namespace Sparrow
} // namespace Scine
