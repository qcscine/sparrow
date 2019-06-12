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
#include <Core/ModuleManager.h>
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/GeometricDerivatives/NormalModeAnalyzer.h>
#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XYZStreamHandler.h>
#include <Utils/IO/FormattedIOUtils.h>
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

  methodWrapper_->setRequiredProperties(requiredProperties);
}

void CalculationHandler::assignSettings() {
  commandLineOptions_.updateSettings(methodWrapper_->settings());
  std::ifstream structureFile(commandLineOptions_.getStructureCoordinatesFile());

  if (structureFile.is_open()) {
    auto structure = Utils::XYZStreamHandler::read(structureFile);
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
}

void CalculationHandler::printResultsToFile() const {
  const std::string& description = results_.getDescription();
  const std::string programCall = commandLineOptions_.getCallStatement();
  const std::string method = commandLineOptions_.getSelectedMethodName();
  std::string conditionalUnderscore = (description == "") ? "" : "_";

  std::string energyName = description + conditionalUnderscore + "energy.dat";

  std::ofstream energyOut(energyName);
  energyOut << std::defaultfloat << std::setprecision(10);
  printHeader(energyOut, "# ");
  energyOut << "# Energy [hartree]:" << std::endl;
  energyOut << results_.getEnergy() << std::endl;
  energyOut.close();

  if (results_.hasGradients()) {
    std::string gradientsName = description + conditionalUnderscore + "gradients.dat";
    std::ofstream gradientsOut(gradientsName);
    gradientsOut << std::scientific;
    printHeader(gradientsOut, "# ");
    gradientsOut << "# Gradients [hartree/bohr]:" << std::endl;
    gradientsOut << results_.getGradients() << std::endl;
    gradientsOut.close();
  }
  if (results_.hasHessian()) {
    std::string hessianName = description + conditionalUnderscore + "hessian.dat";
    std::ofstream hessianOut(hessianName);
    hessianOut << std::scientific;
    printHeader(hessianOut, "# ");
    hessianOut << "# Hessian [hartree/bohr^2]:" << std::endl;
    hessianOut << results_.getHessian();
    hessianOut.close();
  }
  if (results_.hasBondOrders()) {
    std::string bondOrdersName = description + conditionalUnderscore + "bond_orders.dat";
    std::ofstream bondOrdersOut(bondOrdersName);
    bondOrdersOut << std::scientific;
    printHeader(bondOrdersOut, "# ");
    bondOrdersOut << "# Bond orders [dimensionless]:" << std::endl;
    bondOrdersOut << results_.getBondOrders().getMatrix();
    bondOrdersOut.close();
  }
}

void CalculationHandler::printPrettyResults(std::ostream& out) const {
  out << std::string(80, '=') << std::endl;

  out << std::defaultfloat << std::setprecision(10);
  out << "Calculation: " << results_.getDescription() << std::endl;
  out << std::endl;
  out << "Energy [hartree]:" << std::endl;
  out << results_.getEnergy() << std::endl;
  out << std::setprecision(6) << std::endl;
  if (results_.hasGradients()) {
    out << "Gradients [hartree/bohr]:" << std::endl;
    Utils::matrixPrettyPrint(out, results_.getGradients());
    out << std::endl;
  }
  if (results_.hasHessian()) {
    if (!commandLineOptions_.suppressNormalModes()) {
      printFrequencyAnalysis(out, results_.getHessian());
    }
    else {
      out << "Hessian [hartree/bohr^2]:" << std::endl;
      Utils::matrixPrettyPrint(out, results_.getHessian());
    }
  }
  if (results_.hasBondOrders()) {
    out << "Bond orders [dimensionless]:" << std::endl;
    Utils::matrixPrettyPrint(out, results_.getBondOrders().getMatrix(), 0.05);
    out << std::endl;
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
} // namespace Sparrow
} // namespace Scine
