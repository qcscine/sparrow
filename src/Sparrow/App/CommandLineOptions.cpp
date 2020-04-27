/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "CommandLineOptions.h"
#include <Core/Exceptions.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <boost/program_options.hpp>
#include <iostream>

namespace Scine {
namespace Sparrow {

namespace {
constexpr const char* helpKey = "help";
constexpr const char* outputToFile = "output_to_file";
constexpr const char* structureKey = "structure";
constexpr const char* methodKey = "method";

constexpr const char* calculateGradients = "gradients";
constexpr const char* calculateHessian = "hessian";
constexpr const char* calculateBondOrders = "bond_orders";
constexpr const char* calculationDescription = "description";
constexpr const char* suppressNormalModesOutput = "suppress_normal_modes";
constexpr const char* wavefunctionOutput = "wavefunction";
constexpr const char* thermochemistry = "thermochemistry";
} // namespace

using namespace boost::program_options;

struct CommandLineOptions::Impl {
  boost::program_options::options_description desc_;
  boost::program_options::variables_map vm_;
  std::string callStatement;
};

CommandLineOptions::CommandLineOptions(int argc, char* argv[]) : pImpl_(std::make_unique<Impl>()) {
  pImpl_->callStatement = generateCallStatement(argc, argv);

  auto helpOption = combineNamesForOptions(helpKey, ",h");
  auto outputToFileOption = combineNamesForOptions(outputToFile, ",o");
  auto chargeOption = combineNamesForOptions(Utils::SettingsNames::molecularCharge, ",c");
  auto multiplicityOption = combineNamesForOptions(Utils::SettingsNames::spinMultiplicity, ",s");
  auto unrestrictedOption = combineNamesForOptions(Utils::SettingsNames::unrestrictedCalculation, ",u");
  auto parameterFileOption = combineNamesForOptions(Utils::SettingsNames::parameterFile, ",p");
  auto parameterRootOption = combineNamesForOptions(Utils::SettingsNames::parameterRootDirectory, ",P");
  auto mixerTypeOption = combineNamesForOptions(Utils::SettingsNames::mixer, ",m");
  auto maxIterationsOption = combineNamesForOptions(Utils::SettingsNames::maxIterations, ",i");
  auto scfThresholdOption = combineNamesForOptions(Utils::SettingsNames::selfConsistanceCriterion, ",t");
  auto methodOption = combineNamesForOptions(methodKey, ",M");
  auto structureOption = combineNamesForOptions(structureKey, ",x");
  auto gradientOption = combineNamesForOptions(calculateGradients, ",G");
  auto hessianOption = combineNamesForOptions(calculateHessian, ",H");
  auto thermochemistryOption = combineNamesForOptions(thermochemistry, ",C");
  auto temperatureOption = combineNamesForOptions(Utils::SettingsNames::temperature, ",T");
  auto suppressNormalModesOption = combineNamesForOptions(suppressNormalModesOutput, ",N");
  auto bondOrdersOption = combineNamesForOptions(calculateBondOrders, ",B");
  auto wavefunctionOption = combineNamesForOptions(wavefunctionOutput, ",W");
  auto descriptionOption = combineNamesForOptions(calculationDescription, ",d");
  auto logOption = combineNamesForOptions(Utils::SettingsNames::loggerVerbosity, ",l");

  // clang-format off
  pImpl_->desc_.add_options()
    (helpOption.c_str(), "prints this help message")
    (outputToFileOption.c_str(), "flag setting the saving of the output matrices as files")
    (chargeOption.c_str(), value<int>(), "sets the molecular charge")
    (multiplicityOption.c_str(), value<int>(), "sets the spin multiplicity")
    (unrestrictedOption.c_str(), "runs unrestricted calculation (UHF)")
    (parameterFileOption.c_str(), value<std::string>(), "sets the desired parameter file")
    (parameterRootOption.c_str(), value<std::string>(), "sets the desired resource directory")
    (mixerTypeOption.c_str(), value<std::string>(), "sets the desired SCF mixer for convergence acceleration")
    (maxIterationsOption.c_str(), value<int>(), "sets the maximum number of SCF iterations")
    (scfThresholdOption.c_str(), value<double>(), "sets the SCF convergence threshold (in hartree)")
    (methodOption.c_str(), value<std::string>(), "sets the method used to perform the calculation")
    (structureOption.c_str(), value<std::string>()->value_name("structure"), "sets the path to the coordinates file")
    (gradientOption.c_str(), "activates the calculation of the atomic gradients")
    (hessianOption.c_str(), "activates the calculation of the molecular Hessian matrix")
    (thermochemistryOption.c_str(), "activates the calculation of thermochemical properties.")
    (temperatureOption.c_str(), value<double>(), "defines the temperature for the thermochemical calculation in Kelvin.")
    (suppressNormalModesOption.c_str(), "suppresses the results of the hessian calculation as normal modes analysis")
    (bondOrdersOption.c_str(), "activates the calculation of the bond order matrix")
    (wavefunctionOption.c_str(), "outputs a molden file for the visualization of orbitals")
    (descriptionOption.c_str(), value<std::string>(), "sets a calculation description which will appear in the output")
    (logOption.c_str(), value<std::string>()->default_value("warning"), "sets whether warnings and errors are printed");

  // clang-format on

  // Set options without name to be structures
  positional_options_description p;
  p.add(structureOption.c_str(), -1);

  store(command_line_parser(argc, argv).options(pImpl_->desc_).positional(p).run(), pImpl_->vm_);
  notify(pImpl_->vm_);
}

CommandLineOptions::~CommandLineOptions() = default;

std::string CommandLineOptions::getCallStatement() const {
  return pImpl_->callStatement;
}

bool CommandLineOptions::helpRequired() const {
  return pImpl_->vm_.count(helpKey) > 0;
}

bool CommandLineOptions::outputToFileRequired() const {
  return pImpl_->vm_.count(outputToFile) > 0;
}

std::string CommandLineOptions::getSelectedMethodName() const {
  if (pImpl_->vm_.count(methodKey) > 0) {
    std::string result = pImpl_->vm_[methodKey].as<std::string>();
    for (auto& c : result)
      c = static_cast<char>(std::toupper(c));
    return result;
  }
  return "PM6";
}

std::string CommandLineOptions::getCalculationDescription() const {
  if (pImpl_->vm_.count(calculationDescription) > 0) {
    return pImpl_->vm_[calculationDescription].as<std::string>();
  }
  return "";
}

std::string CommandLineOptions::getLoggerVerbosity() const {
  return pImpl_->vm_[Utils::SettingsNames::loggerVerbosity].as<std::string>();
}

bool CommandLineOptions::gradientRequired() const {
  return pImpl_->vm_.count(calculateGradients) > 0;
}

bool CommandLineOptions::hessianRequired() const {
  return pImpl_->vm_.count(calculateHessian) > 0;
}

bool CommandLineOptions::bondOrdersRequired() const {
  return pImpl_->vm_.count(calculateBondOrders) > 0;
}

bool CommandLineOptions::suppressNormalModes() const {
  return pImpl_->vm_.count(suppressNormalModesOutput) > 0;
}

bool CommandLineOptions::wavefunctionRequired() const {
  return pImpl_->vm_.count(wavefunctionOutput) > 0;
}

bool CommandLineOptions::thermochemistryRequired() const {
  return pImpl_->vm_.count(thermochemistry) > 0;
}

template<class CharPtrType, class StringType>
std::string CommandLineOptions::combineNamesForOptions(CharPtrType nameOfOption, StringType abbreviatedOption) const {
  auto combinedString = static_cast<std::string>(nameOfOption) + static_cast<std::string>(abbreviatedOption);
  return combinedString;
}

void CommandLineOptions::updateSettings(Utils::Settings& settingsToUpdate) const {
  using namespace Utils::SettingsNames;
  if (validOptionToSet(molecularCharge, settingsToUpdate)) {
    settingsToUpdate.modifyInt(molecularCharge, pImpl_->vm_[molecularCharge].as<int>());
  }
  if (validOptionToSet(spinMultiplicity, settingsToUpdate)) {
    settingsToUpdate.modifyInt(spinMultiplicity, pImpl_->vm_[spinMultiplicity].as<int>());
  }
  if (validOptionToSet(unrestrictedCalculation, settingsToUpdate)) {
    settingsToUpdate.modifyBool(unrestrictedCalculation, true);
  }
  if (validOptionToSet(mixer, settingsToUpdate)) {
    settingsToUpdate.modifyString(mixer, pImpl_->vm_[mixer].as<std::string>());
  }
  if (validOptionToSet(maxIterations, settingsToUpdate)) {
    settingsToUpdate.modifyInt(maxIterations, pImpl_->vm_[maxIterations].as<int>());
  }
  if (validOptionToSet(selfConsistanceCriterion, settingsToUpdate)) {
    settingsToUpdate.modifyDouble(selfConsistanceCriterion, pImpl_->vm_[selfConsistanceCriterion].as<double>());
  }
  if (validOptionToSet(parameterFile, settingsToUpdate)) {
    settingsToUpdate.modifyString(parameterFile, pImpl_->vm_[parameterFile].as<std::string>());
  }
  if (validOptionToSet(parameterRootDirectory, settingsToUpdate)) {
    settingsToUpdate.modifyString(parameterRootDirectory, pImpl_->vm_[parameterRootDirectory].as<std::string>());
  }
  if (validOptionToSet(temperature, settingsToUpdate)) {
    settingsToUpdate.modifyDouble(temperature, pImpl_->vm_[temperature].as<double>());
  }
}

template<class CharType>
bool CommandLineOptions::validOptionToSet(CharType optionIdentifier, const Utils::Settings& settings) const {
  return pImpl_->vm_.count(optionIdentifier) && settings.valueExists(optionIdentifier);
}

void CommandLineOptions::printHelp(std::ostream& out) const {
  out << "Usage: <executable> [options]\n";
  out << pImpl_->desc_;
}

std::string CommandLineOptions::getStructureCoordinatesFile() const {
  if (pImpl_->vm_.count(structureKey) > 0) {
    return pImpl_->vm_[structureKey].as<std::string>();
  }
  throw Core::InitializationException("No coordinate file given.");
}

std::string CommandLineOptions::generateCallStatement(int argc, char* argv[]) const {
  std::string callStatement = "";
  for (int i = 0; i < argc; ++i) {
    callStatement += " " + std::string(argv[i]);
  }
  return callStatement;
}

} // namespace Sparrow
} // namespace Scine
