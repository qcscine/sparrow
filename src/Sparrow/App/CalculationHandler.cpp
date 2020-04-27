/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Internal Includes */
#include "CalculationHandler.h"
#include "CommandLineOptions.h"
#include "SparrowInitializer.h"
/* External Includes */
#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/WavefunctionOutputGenerator.h>
#include <Core/ModuleManager.h>
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/GeometricDerivatives/NormalModeAnalyzer.h>
#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/FormattedIOUtils.h>
#include <Utils/IO/Logger.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <boost/asio/ip/host_name.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

namespace Scine {
namespace Sparrow {

CalculationHandler::CalculationHandler(CommandLineOptions& options, SparrowInitializer& initializer)
  : commandLineOptions_(options) {
  Utils::Log::startConsoleLogging(commandLineOptions_.getLoggerVerbosity());
  try {
    methodWrapper_ = initializer.getManager().get<Core::Calculator>(commandLineOptions_.getSelectedMethodName());
  }
  catch (std::exception& e) {
    throw MethodNotAvailableException(commandLineOptions_.getSelectedMethodName());
  }
  methodWrapper_->settings().modifyString(Utils::SettingsNames::parameterRootDirectory,
                                          SparrowInitializer::getResourceDirectory());
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
  if (commandLineOptions_.bondOrdersRequired()) {
    requiredProperties.addProperty(Utils::Property::BondOrderMatrix);
  }
  if (commandLineOptions_.thermochemistryRequired()) {
    requiredProperties.addProperty(Utils::Property::Thermochemistry);
  }

  methodWrapper_->setRequiredProperties(requiredProperties);
}

void CalculationHandler::assignSettings() {
  commandLineOptions_.updateSettings(methodWrapper_->settings());
  std::ifstream structureFile(commandLineOptions_.getStructureCoordinatesFile());

  if (structureFile.is_open()) {
    auto structure = Utils::XyzStreamHandler::read(structureFile);
    methodWrapper_->setStructure(structure);
  }
  else {
    throw FileInaccessibleException(commandLineOptions_.getStructureCoordinatesFile());
  }
}

void CalculationHandler::printSettings(std::ostream& out, const Utils::Settings& settings, std::string commentChar) const {
  auto keys = settings.getKeys();
  out << commentChar << std::setw(25) << std::left << "Calculation Modus: ";
  std::string modus = "Spin-restricted formalism";
  if (settings.valueExists(Utils::SettingsNames::unrestrictedCalculation)) {
    if (settings.getBool(Utils::SettingsNames::unrestrictedCalculation))
      modus = "Spin-unrestricted formalism";
  }
  out << modus << std::endl;

  for (const auto& key : keys) {
    if (key == Utils::SettingsNames::molecularCharge) {
      out << commentChar << std::setw(25) << std::left << "Molecular charge: ";
      out << settings.getInt(key) << std::endl;
    }
    else if (key == Utils::SettingsNames::spinMultiplicity) {
      out << commentChar << std::setw(25) << std::left << "Spin multiplicity: ";
      out << settings.getInt(key) << std::endl;
    }
    else if (key == Utils::SettingsNames::maxIterations) {
      out << commentChar << std::setw(25) << std::left << "Max Iterations: ";
      out << settings.getInt(key) << std::endl;
    }
    else if (key == Utils::SettingsNames::selfConsistanceCriterion) {
      out << commentChar << std::setw(25) << std::left << "Convergence threshold: ";
      out << settings.getDouble(key) << std::left << " hartree" << std::endl;
    }
    else if (key == Utils::SettingsNames::mixer) {
      out << commentChar << std::setw(25) << std::left << "Convergence accelerator: ";
      out << settings.getString(key) << std::endl;
    }
    else if (key == Utils::SettingsNames::parameterFile) {
      std::string root = settings.getString(Utils::SettingsNames::parameterRootDirectory);
      std::string par = settings.getString(Utils::SettingsNames::parameterFile);
      auto parFull = Utils::NativeFilenames::combinePathSegments(root, par);
      par += settings.getString(Utils::SettingsNames::parameterFile);
      out << commentChar << std::setw(25) << std::left << "Parameters: ";
      out << parFull << std::endl;
    }
    else if (key == Utils::SettingsNames::temperature && commandLineOptions_.thermochemistryRequired()) {
      out << commentChar << std::setw(25) << std::left << "Temperature: ";
      out << settings.getDouble(key) << std::endl;
    }
  }
}

void CalculationHandler::calculate(std::ostream& out) {
  printHeader(out);
  start_ = std::chrono::system_clock::now();
  results_ = methodWrapper_->calculate(commandLineOptions_.getCalculationDescription());
  end_ = std::chrono::system_clock::now();
  printPrettyResults(out);
  printFooter(out);

  if (commandLineOptions_.outputToFileRequired())
    printResultsToFile();
  if (commandLineOptions_.wavefunctionRequired()) {
    printWavefunction();
  }
}

void CalculationHandler::printResultsToFile() const {
  const std::string& description = results_.get<Utils::Property::Description>();
  const std::string programCall = commandLineOptions_.getCallStatement();
  const std::string method = commandLineOptions_.getSelectedMethodName();
  std::string conditionalUnderscore = (description == "") ? "" : "_";

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
  if (results_.has<Utils::Property::BondOrderMatrix>()) {
    std::string bondOrdersName = description + conditionalUnderscore + "bond_orders.dat";
    std::ofstream bondOrdersOut(bondOrdersName);
    bondOrdersOut << std::scientific;
    printHeader(bondOrdersOut, "# ");
    bondOrdersOut << "# Bond orders [dimensionless]:" << std::endl;
    bondOrdersOut << results_.get<Utils::Property::BondOrderMatrix>().getMatrix();
    bondOrdersOut.close();
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
  if (results_.has<Utils::Property::BondOrderMatrix>()) {
    out << "Bond orders [dimensionless]:" << std::endl;
    Utils::matrixPrettyPrint(out, results_.get<Utils::Property::BondOrderMatrix>().getMatrix(), 0.05);
    out << std::endl;
  }
  if (results_.has<Utils::Property::Thermochemistry>()) {
    out << "Thermochemistry:" << std::endl;
    Utils::prettyPrint(out, results_.get<Utils::Property::Thermochemistry>() * Utils::Constants::kJPerMol_per_hartree);
  }
  out << std::string(80, '=') << std::endl;
}

void CalculationHandler::printHeader(std::ostream& out, std::string commentChar) const {
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
  out << std::endl;
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_ - start_);
  out << "Total job duration: " << duration.count() << " milliseconds." << std::endl;
  out << std::endl;
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

  Utils::NormalModeAnalyzer frequencyAnalyzer(matrix, elements, positions);
  auto normalModesContainer = frequencyAnalyzer.calculateNormalModes();
  Utils::matrixPrettyPrint(out, normalModesContainer, elements);
}

void CalculationHandler::printWavefunction() const {
  using namespace Utils::SettingsNames;
  auto& manager = Core::ModuleManager::getInstance();
  auto wfGenerator = manager.get<Core::WavefunctionOutputGenerator>(commandLineOptions_.getSelectedMethodName());

  if (!methodWrapper_->results().has<Utils::Property::SuccessfulCalculation>() ||
      !methodWrapper_->results().get<Utils::Property::SuccessfulCalculation>())
    throw std::runtime_error("Wavefunction requested, but ground state calculation was not successful.");

  wfGenerator->settings() = methodWrapper_->settings();
  wfGenerator->setStructure(*methodWrapper_->getStructure());
  wfGenerator->loadState(methodWrapper_->getState());

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
