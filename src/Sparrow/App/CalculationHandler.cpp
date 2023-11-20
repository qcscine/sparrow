/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Internal Includes */
#include "CalculationHandler.h"
#include "CommandLineOptions.h"
/* External Includes */
#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/WavefunctionOutputGenerator.h>
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/GeometricDerivatives/NormalModeAnalysis.h>
#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/FormattedIOUtils.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/Math/AutomaticDifferentiation/Second3D.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <boost/asio/ip/host_name.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

namespace Scine {
namespace Sparrow {

CalculationHandler::CalculationHandler(CommandLineOptions& options) : commandLineOptions_(options) {
  try {
    methodWrapper_ = Core::ModuleManager::getInstance().get<Core::Calculator>(commandLineOptions_.getSelectedMethodName(),
                                                                              "Sparrow");
  }
  catch (std::exception& e) {
    throw MethodNotAvailableException(commandLineOptions_.getSelectedMethodName());
  }

  if (commandLineOptions_.excitedStatesRequired()) {
    std::string excitedStatesMethod;
    if (isNDDO()) {
      excitedStatesMethod = "CIS-NDDO";
    }
    else if (isDFTB()) {
      excitedStatesMethod = "TD-DFTB";
    }
    else {
      throw MethodNotAvailableException(commandLineOptions_.getSelectedMethodName());
    }

    try {
      excitedStatesCalculator_ = Core::ModuleManager::getInstance().get<Core::CalculatorWithReference>(excitedStatesMethod);
    }
    catch (...) {
      throw MethodNotAvailableException(excitedStatesMethod);
    }
  }

  if (commandLineOptions_.orbitalSteeringRequired()) {
    try {
      orbitalSteerer_ = Core::ModuleManager::getInstance().get<Core::CalculatorWithReference>("orbital_steering");
      orbitalSteerer_->setReferenceCalculator(methodWrapper_);
    }
    catch (...) {
      throw MethodNotAvailableException("orbital_steering");
    }
  }
  assignPropertiesToCalculate();
  assignSettings();
}

void CalculationHandler::assignPropertiesToCalculate() {
  Utils::PropertyList requiredProperties;
  if (commandLineOptions_.gradientRequired()) {
    requiredProperties.addProperty(Utils::Property::Gradients);
  }
  if (commandLineOptions_.hessianRequired()) {
    requiredProperties.addProperty(Utils::Property::Gradients);
    requiredProperties.addProperty(Utils::Property::Hessian);
  }
  if (commandLineOptions_.atomicHessiansRequired()) {
    requiredProperties.addProperty(Utils::Property::AtomicHessians);
  }
  if (commandLineOptions_.bondOrdersRequired()) {
    requiredProperties.addProperty(Utils::Property::BondOrderMatrix);
  }
  if (commandLineOptions_.excitedStatesRequired()) {
    requiredProperties.addProperty(Utils::Property::DipoleMatrixMO);
  }
  if (commandLineOptions_.thermochemistryRequired()) {
    requiredProperties.addProperty(Utils::Property::Thermochemistry);
  }

  methodWrapper_->setRequiredProperties(requiredProperties);
}

void CalculationHandler::assignSettings() {
  commandLineOptions_.updateSettings(methodWrapper_->settings());
  commandLineOptions_.updateLogger(methodWrapper_->getLog());
  std::ifstream structureFile(commandLineOptions_.getStructureCoordinatesFile());

  if (structureFile.is_open()) {
    auto structure = Utils::XyzStreamHandler::read(structureFile);
    methodWrapper_->setStructure(structure);
  }
  else {
    throw FileInaccessibleException(commandLineOptions_.getStructureCoordinatesFile());
  }
  if (excitedStatesCalculator_) {
    commandLineOptions_.updateLogger(excitedStatesCalculator_->getLog());
    commandLineOptions_.updateExcitedStatesSettings(excitedStatesCalculator_->settings());
    excitedStatesCalculator_->applySettings();
  }
  if (orbitalSteerer_) {
    commandLineOptions_.updateLogger(orbitalSteerer_->getLog());
    commandLineOptions_.updateOrbitalSteeringSettings(orbitalSteerer_->settings());
    orbitalSteerer_->applySettings();
  }
}

void CalculationHandler::printSettings(std::ostream& out, const Utils::Settings& settings, const std::string& commentChar) const {
  auto keys = settings.getKeys();
  out << commentChar << std::setw(25) << std::left << "Calculation Modus: ";
  std::string modus = "Spin-restricted formalism";
  if (settings.valueExists(Utils::SettingsNames::spinMode)) {
    if (settings.getString(Utils::SettingsNames::spinMode) == "unrestricted") {
      modus = "Spin-unrestricted formalism";
    }
  }
  out << modus << std::endl;

  for (const auto& key : keys) {
    if (key == Utils::SettingsNames::molecularCharge) {
      out << commentChar << std::setw(27) << std::left << "Molecular charge: ";
      out << settings.getInt(key) << std::endl;
    }
    else if (key == Utils::SettingsNames::spinMultiplicity) {
      out << commentChar << std::setw(27) << std::left << "Spin multiplicity: ";
      out << settings.getInt(key) << std::endl;
    }
    else if (key == Utils::SettingsNames::maxScfIterations) {
      out << commentChar << std::setw(27) << std::left << "Max Iterations: ";
      out << settings.getInt(key) << std::endl;
    }
    else if (key == Utils::SettingsNames::selfConsistenceCriterion) {
      out << commentChar << std::setw(27) << std::left << "Convergence threshold: ";
      out << settings.getDouble(key) << std::left << " hartree" << std::endl;
    }
    else if (key == Utils::SettingsNames::mixer) {
      out << commentChar << std::setw(27) << std::left << "Convergence accelerator: ";
      out << settings.getString(key) << std::endl;
    }
    else if (key == Utils::SettingsNames::methodParameters) {
      const std::string par = settings.getString(Utils::SettingsNames::methodParameters);
      out << commentChar << std::setw(27) << std::left << "Parameters: ";
      out << par << std::endl;
    }
    else if (key == Utils::SettingsNames::temperature && commandLineOptions_.thermochemistryRequired()) {
      out << commentChar << std::setw(27) << std::left << "Temperature: ";
      out << settings.getDouble(key) << std::endl;
    }
    else if (key == Utils::SettingsNames::symmetryNumber && commandLineOptions_.thermochemistryRequired()) {
      out << commentChar << std::setw(27) << std::left << "Molecular symmetry number: ";
      out << settings.getInt(key) << std::endl;
    }
  }
  out << std::endl;

  if (excitedStatesCalculator_) {
    out << commentChar << std::setw(25) << std::left << "Excited states: ";
    out << excitedStatesCalculator_->name() << std::endl;
    out << commentChar << std::setw(25) << "Number of roots: ";
    out << excitedStatesCalculator_->settings().getInt(Utils::SettingsNames::numberOfEigenstates) << std::endl;
    out << commentChar << std::setw(25) << "Initial guess space: ";
    out << excitedStatesCalculator_->settings().getInt(Utils::SettingsNames::initialSubspaceDimension) << std::endl;
  }
  if (orbitalSteerer_) {
    out << commentChar << std::setw(40) << std::left << "Orbital Steerer: " << std::endl;
    out << commentChar << std::setw(40) << "Max iterations for steerer: ";
    out << orbitalSteerer_->settings().getInt(Utils::SettingsNames::maxScfIterations) << std::endl;
    out << commentChar << std::setw(40) << "SCF mixer for steerer: ";
    out << orbitalSteerer_->settings().getString(Utils::SettingsNames::mixer) << std::endl;
    out << commentChar << std::setw(40) << "Number of orbitals to mix: ";
    out << orbitalSteerer_->settings().getInt("number_orbitals_to_mix") << std::endl;
    out << commentChar << std::setw(40) << "Number of orbitals to consider: ";
    int nOrbsToConsider = orbitalSteerer_->settings().getInt("number_orbitals_to_consider");
    out << (nOrbsToConsider == 0 ? "All" : std::to_string(nOrbsToConsider)) << std::endl;
    out << commentChar << std::setw(40) << "Minimal mixing angle: ";
    out << orbitalSteerer_->settings().getDouble("minimal_mixing_angle") << std::endl;
    out << commentChar << std::setw(40) << "Maximal mixing angle: ";
    out << orbitalSteerer_->settings().getDouble("maximal_mixing_angle") << std::endl;
    out << commentChar << std::setw(40) << "Number of repeated mixes: ";
    out << commandLineOptions_.getNumberOfOrbitalSteers() << std::endl;
  }
}

void CalculationHandler::calculate(std::ostream& out) {
  printHeader(out);
  auto startGS = std::chrono::system_clock::now();
  results_ = methodWrapper_->calculate(commandLineOptions_.getCalculationDescription());
  auto endGS = std::chrono::system_clock::now();
  printCalculationConverged(out);
  groundStateTime_ = std::chrono::duration_cast<std::chrono::milliseconds>(endGS - startGS).count();

  if (commandLineOptions_.orbitalSteeringRequired()) {
    auto startES = std::chrono::system_clock::now();
    for (int i = 0; i < commandLineOptions_.getNumberOfOrbitalSteers(); ++i) {
      results_ = orbitalSteerer_->calculate();
    }
    auto endES = std::chrono::system_clock::now();
    excitedStateTime_ = std::chrono::duration_cast<std::chrono::milliseconds>(endES - startES).count();
  }

  if (commandLineOptions_.excitedStatesRequired()) {
    excitedStatesCalculator_->setReferenceCalculator(methodWrapper_);
    auto startES = std::chrono::system_clock::now();
    excitedStatesResults_ = excitedStatesCalculator_->calculate();
    auto endES = std::chrono::system_clock::now();
    excitedStateTime_ = std::chrono::duration_cast<std::chrono::milliseconds>(endES - startES).count();
  }

  printPrettyResults(out);
  printFooter(out);

  if (commandLineOptions_.outputToFileRequired())
    printResultsToFile();
  if (commandLineOptions_.wavefunctionRequired()) {
    printWavefunction();
  }
}

void CalculationHandler::printCalculationConverged(std::ostream& out) {
  if (results_.get<Utils::Property::SuccessfulCalculation>()) {
    out << "SCF converged!" << std::endl;
  }
  else {
    out << "SCF not converged!" << std::endl;
  }
}

void CalculationHandler::printResultsToFile() const {
  const std::string& description = results_.get<Utils::Property::Description>();
  const std::string programCall = commandLineOptions_.getCallStatement();
  const std::string method = commandLineOptions_.getSelectedMethodName();
  std::string conditionalUnderscore = (description.empty()) ? "" : "_";

  std::string energyName = description + conditionalUnderscore + "energy.dat";

  std::ofstream energyOut(energyName);
  energyOut << std::defaultfloat << std::setprecision(10);
  printHeader(energyOut, "# ");
  energyOut << "# Energy [hartree]:" << std::endl;
  energyOut << results_.get<Utils::Property::Energy>() << std::endl;
  energyOut.close();

  if (results_.has<Utils::Property::Gradients>()) {
    std::string gradientsName = description + conditionalUnderscore + "gradients.dat";
    std::ofstream gradientsOut(gradientsName);
    gradientsOut << std::scientific;
    printHeader(gradientsOut, "# ");
    gradientsOut << "# Gradients [hartree/bohr]:" << std::endl;
    gradientsOut << results_.get<Utils::Property::Gradients>() << std::endl;
    gradientsOut.close();
  }
  if (results_.has<Utils::Property::Hessian>()) {
    std::string hessianName = description + conditionalUnderscore + "hessian.dat";
    std::ofstream hessianOut(hessianName);
    hessianOut << std::scientific;
    printHeader(hessianOut, "# ");
    hessianOut << "# Hessian [hartree/bohr^2]:" << std::endl;
    hessianOut << results_.get<Utils::Property::Hessian>();
    hessianOut.close();
  }
  if (results_.has<Utils::Property::AtomicHessians>()) {
    std::string atomicHessiansName = description + conditionalUnderscore + "atomic_hessians.dat";
    std::ofstream atomicHessiansOut(atomicHessiansName);
    atomicHessiansOut << std::scientific;
    printHeader(atomicHessiansOut, "# ");
    atomicHessiansOut << "# Atomic Hessian [hartree/bohr^2]:" << std::endl;

    int i = 1;
    for (Utils::AutomaticDifferentiation::Second3D ah : results_.get<Utils::Property::AtomicHessians>()) {
      atomicHessiansOut << std::right << std::scientific << std::setprecision(10);
      atomicHessiansOut << std::setw(15) << "Atomic Index:" << std::setw(5) << i << std::endl;
      atomicHessiansOut << std::setw(30) << "X" << std::setw(20) << "Y" << std::setw(20) << "Z" << std::endl;
      atomicHessiansOut << std::setw(10) << "X " << std::setw(20) << ah.XX() << std::setw(20) << ah.XY()
                        << std::setw(20) << ah.XZ() << std::endl;
      atomicHessiansOut << std::setw(10) << "Y " << std::setw(20) << ah.YX() << std::setw(20) << ah.YY()
                        << std::setw(20) << ah.YZ() << std::endl;
      atomicHessiansOut << std::setw(10) << "Z " << std::setw(20) << ah.ZX() << std::setw(20) << ah.ZY()
                        << std::setw(20) << ah.ZZ() << std::endl;
      atomicHessiansOut << std::endl;
      i++;
    }
    atomicHessiansOut.close();
  }
  if (results_.has<Utils::Property::BondOrderMatrix>()) {
    std::string bondOrdersName = description + conditionalUnderscore + "bond_orders.dat";
    std::ofstream bondOrdersOut(bondOrdersName);
    bondOrdersOut << std::scientific;
    printHeader(bondOrdersOut, "# ");
    bondOrdersOut << "# Bond orders [dimensionless]:" << std::endl;
    bondOrdersOut << results_.get<Utils::Property::BondOrderMatrix>().getMatrix();
    bondOrdersOut.close();
  }
  if (excitedStatesResults_.has<Utils::Property::ExcitedStates>()) {
    enum class CISCalculationType { unrestricted, singlet, triplet };
    auto printLambda = [&](CISCalculationType type) {
      const auto& excitedStates = type == CISCalculationType::singlet
                                      ? excitedStatesResults_.get<Utils::Property::ExcitedStates>().singlet
                                      : type == CISCalculationType::triplet
                                            ? excitedStatesResults_.get<Utils::Property::ExcitedStates>().triplet
                                            : excitedStatesResults_.get<Utils::Property::ExcitedStates>().unrestricted;

      std::string multiplicity = type == CISCalculationType::singlet ? "singlet" : "triplet";
      if (type == CISCalculationType::unrestricted)
        multiplicity = "unrestricted";
      std::string transitionEnergiesName = description + conditionalUnderscore + excitedStatesCalculator_->name() +
                                           "_" + multiplicity + "_vertical_transition_energy.dat";
      std::string transitionDipolesName = description + conditionalUnderscore + excitedStatesCalculator_->name() + "_" +
                                          multiplicity + "_transition_dipoles.dat";
      std::string oscillatorStrengthName = description + conditionalUnderscore + excitedStatesCalculator_->name() +
                                           "_" + multiplicity + "_oscillator_strength.dat";
      std::string excitedStatesName = description + conditionalUnderscore + excitedStatesCalculator_->name() + "_" +
                                      multiplicity + "_excited_states.dat";
      std::ofstream verticalTransitionEnergyOut(transitionEnergiesName);
      std::ofstream transitionDipolesOut(transitionDipolesName);
      std::ofstream oscillatorStrengthOut(oscillatorStrengthName);
      std::ofstream excitedStatesOut(excitedStatesName);
      excitedStatesOut << std::scientific;
      transitionDipolesOut << std::scientific;
      oscillatorStrengthOut << std::scientific;
      verticalTransitionEnergyOut << std::scientific;
      printHeader(excitedStatesOut, "# ");
      printHeader(transitionDipolesOut, "# ");
      printHeader(oscillatorStrengthOut, "# ");
      printHeader(verticalTransitionEnergyOut, "# ");
      verticalTransitionEnergyOut << "# " + multiplicity + " vertical transition energies [hartree]" << std::endl;
      verticalTransitionEnergyOut << excitedStates->eigenStates.eigenValues << std::endl;
      transitionDipolesOut << "# " + multiplicity + " transition dipoles [a.u.]" << std::endl;
      transitionDipolesOut << excitedStates->transitionDipoles << std::endl;
      oscillatorStrengthOut << "# " + multiplicity + " oscillator strengths " << std::endl;
      oscillatorStrengthOut << Utils::TransitionDipoleCalculator::transitionDipoleMomentToOscillatorStrength(
                                   excitedStates->transitionDipoles, excitedStates->eigenStates.eigenValues)
                            << std::endl;
      excitedStatesOut << "# " + multiplicity + " excited states coefficients" << std::endl;
      excitedStatesOut << excitedStates->eigenStates.eigenVectors << std::endl;
    };
    if (excitedStatesResults_.get<Utils::Property::ExcitedStates>().singlet) {
      printLambda(CISCalculationType::singlet);
    }
    if (excitedStatesResults_.get<Utils::Property::ExcitedStates>().triplet) {
      printLambda(CISCalculationType::triplet);
    }
    if (excitedStatesResults_.get<Utils::Property::ExcitedStates>().unrestricted) {
      printLambda(CISCalculationType::unrestricted);
    }
  }
  if (results_.has<Utils::Property::Thermochemistry>()) {
    std::string thermoName = description + conditionalUnderscore + "thermochemistry.dat";
    std::ofstream thermoOut(thermoName);
    thermoOut << std::scientific;
    printHeader(thermoOut, "# ");
    thermoOut << "# Thermochemistry:" << std::endl;
    Utils::prettyPrint(thermoOut, results_.get<Utils::Property::Thermochemistry>() * Utils::Constants::kJPerMol_per_hartree);
    thermoOut.close();
  }
}

void CalculationHandler::printPrettyResults(std::ostream& out) const {
  out << std::string(80, '=') << std::endl;
  if (!results_.get<Utils::Property::SuccessfulCalculation>()) {
    out << "Be warned: the scf method has not converged." << std::endl;
  }

  out << std::defaultfloat << std::setprecision(10);
  out << "Calculation: " << results_.get<Utils::Property::Description>() << std::endl;
  out << std::endl;
  out << "Energy [hartree]:" << std::endl;
  out << results_.get<Utils::Property::Energy>() << std::endl;
  out << std::setprecision(6) << std::endl;

  if (results_.has<Utils::Property::Gradients>()) {
    out << "Gradients [hartree/bohr]:" << std::endl;
    Utils::matrixPrettyPrint(out, results_.get<Utils::Property::Gradients>());
    out << std::endl;
  }

  if (results_.has<Utils::Property::Hessian>()) {
    if (!commandLineOptions_.suppressNormalModes()) {
      printFrequencyAnalysis(out, results_.get<Utils::Property::Hessian>());
    }
    else {
      out << "Hessian [hartree/bohr^2]:" << std::endl;
      Utils::matrixPrettyPrint(out, results_.get<Utils::Property::Hessian>());
    }
  }
  if (results_.has<Utils::Property::AtomicHessians>()) {
    out << "Atomic Hessians [hartree/bohr^2]:" << std::endl;
    int i = 1;
    for (Utils::AutomaticDifferentiation::Second3D ah : results_.get<Utils::Property::AtomicHessians>()) {
      out << std::right << std::scientific << std::setprecision(10);
      out << std::setw(15) << "Atomic Index:" << std::setw(5) << i << std::endl;
      out << std::setw(30) << "X" << std::setw(20) << "Y" << std::setw(20) << "Z" << std::endl;
      out << std::setw(10) << "X " << std::setw(20) << ah.XX() << std::setw(20) << ah.XY() << std::setw(20) << ah.XZ()
          << std::endl;
      out << std::setw(10) << "Y " << std::setw(20) << ah.YX() << std::setw(20) << ah.YY() << std::setw(20) << ah.YZ()
          << std::endl;
      out << std::setw(10) << "Z " << std::setw(20) << ah.ZX() << std::setw(20) << ah.ZY() << std::setw(20) << ah.ZZ()
          << std::endl;
      out << std::endl;
      i++;
    }
  }
  if (results_.has<Utils::Property::BondOrderMatrix>()) {
    out << "Bond orders [dimensionless]:" << std::endl;
    Utils::matrixPrettyPrint(out, results_.get<Utils::Property::BondOrderMatrix>().getMatrix(), 0.05);
    out << std::endl;
  }

  if (excitedStatesResults_.has<Utils::Property::ExcitedStates>()) {
    printExcitedStates(out, excitedStatesResults_.get<Utils::Property::ExcitedStates>());
  }

  if (results_.has<Utils::Property::Thermochemistry>()) {
    out << "Thermochemistry:" << std::endl;
    Utils::prettyPrint(out, results_.get<Utils::Property::Thermochemistry>() * Utils::Constants::kJPerMol_per_hartree);
  }
  out << std::string(80, '=') << std::endl;
}

void CalculationHandler::printHeader(std::ostream& out, const std::string& commentChar) const {
  out << commentChar + "SPARROW, command-line program\n";
  out << commentChar + "Host: " << boost::asio::ip::host_name() << std::endl;
  out << commentChar + "Start: ";
  printTime(out);
  out << commentChar << std::endl;
  out << commentChar + "Program call: " + commandLineOptions_.getCallStatement() << std::endl;
  out << commentChar << std::setw(25) << std::left << "Method: ";
  out << commandLineOptions_.getSelectedMethodName() << std::endl;
  printSettings(out, methodWrapper_->settings(), commentChar);
  out << commentChar << std::endl;
}

void CalculationHandler::printFooter(std::ostream& out) const {
  out << std::fixed << std::setprecision(0) << std::endl;
  out << std::left << std::setw(40) << "Ground state calculation time: " << std::right << std::setw(5)
      << groundStateTime_ << " ms." << std::endl;
  if (results_.has<Utils::Property::Hessian>() && (!commandLineOptions_.suppressNormalModes())) {
    out << std::left << std::setw(40)
        << "Frequency analysis (Hessian diagonalization + normal modes construction) time: " << std::right
        << std::setw(5) << hessianDiagTime_ << " ms." << std::endl;
  }
  if (excitedStatesResults_.has<Utils::Property::ExcitedStates>()) {
    out << std::left << std::setw(40) << "Excited states calculation time: " << std::right << std::setw(5)
        << excitedStateTime_ << " ms." << std::endl;
  }
  if (orbitalSteerer_) {
    out << std::left << std::setw(40) << "Orbital Steerer calculation time: " << std::right << std::setw(5)
        << excitedStateTime_ << " ms." << std::endl;
  }
  out << std::left << std::setw(40) << "Total job duration: " << std::right << std::setw(5)
      << groundStateTime_ + excitedStateTime_ + hessianDiagTime_ << " ms." << std::endl;
  out << std::left << std::endl;
  out << "End: ";
  printTime(out);
}

void CalculationHandler::printTime(std::ostream& out) const {
  auto timePoint = std::chrono::system_clock::now();
  auto time_t = std::chrono::system_clock::to_time_t(timePoint);
  out << ctime(&time_t);
}

void CalculationHandler::printFrequencyAnalysis(std::ostream& out, const Utils::HessianMatrix& matrix) const {
  out << "Frequency Analysis: frequencies [cm^-1] and mass-weighted normal modes." << std::endl << std::endl;
  const auto elements = methodWrapper_->getStructure()->getElements();
  const auto& positions = methodWrapper_->getPositions();

  auto start = std::chrono::system_clock::now();
  auto normalModesContainer = Utils::NormalModeAnalysis::calculateNormalModes(matrix, elements, positions);
  auto end = std::chrono::system_clock::now();

  Utils::matrixPrettyPrint(out, normalModesContainer, elements);
  hessianDiagTime_ = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
}

void CalculationHandler::printExcitedStates(std::ostream& out, const Utils::SpinAdaptedElectronicTransitionResult& matrix) const {
  std::string method = isNDDO() ? "NDDO-CIS" : "TD-DFTB";
  out << "Excited states calculation with " + method + ".\n" << std::endl;
  Utils::prettyPrint(out, matrix);
}

bool CalculationHandler::isNDDO() const {
  return (commandLineOptions_.getSelectedMethodName() == "PM6" ||
          commandLineOptions_.getSelectedMethodName() == "AM1" || commandLineOptions_.getSelectedMethodName() == "RM1" ||
          commandLineOptions_.getSelectedMethodName() == "PM3" || commandLineOptions_.getSelectedMethodName() == "MNDO");
}

bool CalculationHandler::isDFTB() const {
  return (commandLineOptions_.getSelectedMethodName() == "DFTB0" || commandLineOptions_.getSelectedMethodName() == "DFTB2" ||
          commandLineOptions_.getSelectedMethodName() == "DFTB3");
}

void CalculationHandler::printWavefunction() const {
  using namespace Utils::SettingsNames;
  auto wfGenerator = std::dynamic_pointer_cast<Core::WavefunctionOutputGenerator>(methodWrapper_);

  if (!methodWrapper_->results().has<Utils::Property::SuccessfulCalculation>() ||
      !methodWrapper_->results().get<Utils::Property::SuccessfulCalculation>())
    throw std::runtime_error("Wavefunction requested, but ground state calculation was not successful.");

  const std::string& description = results_.get<Utils::Property::Description>();
  std::string conditionalUnderscore = description.empty() ? "" : "_";

  std::string filename = description + conditionalUnderscore + "wavefunction.molden.input";
  std::ofstream moldenOut(filename);
  if (!moldenOut.is_open()) {
    throw std::runtime_error("Error opening file " + filename);
  }
  wfGenerator->generateWavefunctionInformation(moldenOut);
}

} // namespace Sparrow
} // namespace Scine
