/**
 * @file Calculator.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "OrbitalSteeringCalculator.h"
#include "OrbitalSteeringSettings.h"
#include "Utils/CalculatorBasics/PropertyList.h"
#include <Sparrow/Implementations/GenericMethodWrapper.h>
#include <Utils/Constants.h>
#include <Utils/Scf/MethodInterfaces/LcaoMethod.h>
#include <Utils/Scf/OrbitalPerturbation/RandomOrbitalMixer.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <iomanip>

namespace Scine {
namespace Sparrow {

OrbitalSteeringCalculator::InvalidCalculatorTypeForOrbitalSteerer::InvalidCalculatorTypeForOrbitalSteerer() {
  error = "Calculator is not a Sparrow method and cannot support Orbital Steering.";
}

OrbitalSteeringCalculator::OrbitalSteeringCalculator()
  : settings_(std::make_unique<OrbitalSteeringSettings>()), numberOfCalculations_(0) {
}

void OrbitalSteeringCalculator::setReferenceCalculator(std::shared_ptr<Core::Calculator> referenceCalculator) {
  method_ = std::dynamic_pointer_cast<GenericMethodWrapper>(referenceCalculator);
  if (!method_) {
    throw InvalidCalculatorTypeForOrbitalSteerer();
  }
  applySettings();
}

void OrbitalSteeringCalculator::referenceCalculation() {
  method_->calculate("");
}

Core::Calculator& OrbitalSteeringCalculator::getReferenceCalculator() {
  return *method_;
}

const Core::Calculator& OrbitalSteeringCalculator::getReferenceCalculator() const {
  return *method_;
}

const Utils::Results& OrbitalSteeringCalculator::calculate() {
  applySettings();

  int freq = settings_->getInt(mixingFrequencyKey);

  if (numberOfCalculations_ % freq != 0) {
    referenceCalculation();
  }
  else {
    auto requiredProperties = method_->getRequiredProperties();
    // If method_ has a required Hessian, you do not want to have to calculate a Hessian JUST to
    // optimize the orbitals. Better to just do a ref calculation afterwards.
    method_->setRequiredProperties(Utils::Property::Energy);
    referenceCalculation();

    auto state = method_->getState();

    double energyBeforeSteering = method_->results().get<Utils::Property::Energy>();

    auto scfMixer = method_->settings().getString(Utils::SettingsNames::mixer);
    int maxIterations = method_->settings().getInt(Utils::SettingsNames::maxScfIterations);
    if (settings_->getString(Utils::SettingsNames::mixer) != OrbitalSteeringSettings::SameMixer) {
      method_->settings().modifyString(Utils::SettingsNames::mixer, settings_->getString(Utils::SettingsNames::mixer));
    }
    if (settings_->getInt(Utils::SettingsNames::maxScfIterations) != OrbitalSteeringSettings::SameNumberOfIterations) {
      method_->settings().modifyInt(Utils::SettingsNames::maxScfIterations,
                                    settings_->getInt(Utils::SettingsNames::maxScfIterations));
    }

    auto& molecularOrbitals = method_->getLcaoMethod().getMolecularOrbitals();
    const auto& occ = method_->getLcaoMethod().getElectronicOccupation();
    auto mixer = occ.isRestricted()
                     ? Utils::OrbitalPerturbation::RandomOrbitalMixer(molecularOrbitals, occ.numberRestrictedElectrons())
                     : Utils::OrbitalPerturbation::RandomOrbitalMixer(molecularOrbitals, occ.numberAlphaElectrons(),
                                                                      occ.numberBetaElectrons());

    mixer.setNumberMixes(settings_->getInt(numberOrbitalsToMixKey));
    mixer.setMaximalMixAngle(Utils::Constants::rad_per_degree * settings_->getDouble(maximalAngleKey));
    mixer.setMinimalMixAngle(Utils::Constants::rad_per_degree * settings_->getDouble(minimalAngleKey));

    int orbitalsToConsider = settings_->getInt(numberOrbitalsToConsiderKey);
    if (orbitalsToConsider == OrbitalSteeringSettings::AllOrbitals) {
      mixer.considerAllOrbitals();
    }
    else {
      mixer.setNumberMixes(orbitalsToConsider);
    }

    auto log = numberOfCalculations_ == 0 ? getLog() : Core::Log::silent();
    mixer.mix(log);
    method_->getLcaoMethod().calculateDensity();

    const auto& res = method_->calculate("");
    method_->settings().modifyString(Utils::SettingsNames::mixer, scfMixer);
    method_->settings().modifyInt(Utils::SettingsNames::maxScfIterations, maxIterations);

    double energyAfterSteering = res.get<Utils::Property::Energy>();
    getLog().debug << "Energy Difference: " << energyAfterSteering - energyBeforeSteering << " Ha" << Core::Log::endl;

    if (energyAfterSteering < energyBeforeSteering && res.get<Utils::Property::SuccessfulCalculation>()) {
      logSteering(energyBeforeSteering, energyAfterSteering);
    }
    else {
      method_->loadState(std::move(state));
    }
    method_->setRequiredProperties(requiredProperties);
    referenceCalculation();
  }
  ++numberOfCalculations_;
  return method_->results();
}

std::string OrbitalSteeringCalculator::name() const {
  return model;
}

Utils::Settings& OrbitalSteeringCalculator::settings() {
  return *settings_;
}

const Utils::Settings& OrbitalSteeringCalculator::settings() const {
  return *settings_;
}

void OrbitalSteeringCalculator::applySettings() {
  if (!method_) {
    throw InvalidCalculatorTypeForOrbitalSteerer();
  }
}

Utils::Results& OrbitalSteeringCalculator::results() {
  return method_->results();
}

const Utils::Results& OrbitalSteeringCalculator::results() const {
  return method_->results();
}

void OrbitalSteeringCalculator::logSteering(double oldEnergy, double newEnergy) {
  getLog().output << Core::Log::endl;
  getLog().output << std::setw(1) << "" << std::string(68, '=') << Core::Log::nl;
  getLog().output << std::left << std::setw(12) << ""
                  << "Orbitals Were Steered, New Density Injected." << std::right << Core::Log::endl;
  getLog().output << std::scientific << Core::Log::endl;
  getLog().output << std::setw(1) << "" << std::string(68, '=') << Core::Log::nl;
  getLog().output << std::setw(2) << "#" << std::setw(65) << "" << std::setw(2) << "#" << Core::Log::nl;
  getLog().output << std::setw(2) << "#" << std::setw(20) << "Old Energy [Ha]" << std::setw(20) << "New Energy [Ha]";
  getLog().output << std::setw(20) << "Difference [Ha]" << std::setw(7) << "#" << Core::Log::nl;
  getLog().output << std::setw(2) << "#" << std::setw(20) << oldEnergy << std::setw(20) << newEnergy;
  getLog().output << std::setw(20) << newEnergy - oldEnergy << std::setw(7) << "#" << Core::Log::nl;
  getLog().output << std::setw(2) << "#" << std::setw(65) << "" << std::setw(2) << "#" << Core::Log::nl;
  getLog().output << std::setw(1) << "" << std::string(68, '=') << Core::Log::endl;
}

OrbitalSteeringCalculator::~OrbitalSteeringCalculator() = default;

} // namespace Sparrow
} // namespace Scine
