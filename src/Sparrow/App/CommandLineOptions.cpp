/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "CommandLineOptions.h"
#include <Core/Exceptions.h>
#include <Core/Log.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
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
constexpr const char* logFilename = "log_filename";

constexpr const char* calculateGradients = "gradients";
constexpr const char* calculateHessian = "hessian";
constexpr const char* calculateAtomicHessians = "atomic_hessians";
constexpr const char* calculateBondOrders = "bond_orders";
constexpr const char* calculationDescription = "description";
constexpr const char* suppressNormalModesOutput = "suppress_normal_modes";
constexpr const char* calculateExcitedStates = "excited_states";
constexpr const char* numberOrbitalSteers = "number_orbital_mixes";
constexpr const char* numberOrbitalsToMix = "number_orbitals_to_mix";
constexpr const char* numberOrbitalsToConsider = "number_orbitals_to_consider";
constexpr const char* maximalMixingAngle = "maximal_mixing_angle";
constexpr const char* minimalMixingAngle = "minimal_mixing_angle";
constexpr const char* maxIterationsSteering = "max_iterations_in_steering";
constexpr const char* scfMixerSteering = "scf_mixer_steering";
constexpr const char* distanceThreshold = "distance_threshold";
constexpr const char* wavefunctionOutput = "wavefunction";
constexpr const char* thermochemistry = "thermochemistry";
constexpr const char* unrestrictedCalculation = "unrestricted_calculation";

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
  auto unrestrictedOption = combineNamesForOptions(unrestrictedCalculation, ",u");
  auto parameterFileOption = combineNamesForOptions(Utils::SettingsNames::methodParameters, ",p");
  auto excitedStatesParameterFileOption = combineNamesForOptions(Utils::SettingsNames::excitedStatesParamFile, ",j");
  auto mixerTypeOption = combineNamesForOptions(Utils::SettingsNames::mixer, ",m");
  auto maxScfIterationsOption = combineNamesForOptions(Utils::SettingsNames::maxScfIterations, ",I");
  auto scfThresholdOption = combineNamesForOptions(Utils::SettingsNames::selfConsistenceCriterion, ",t");
  auto methodOption = combineNamesForOptions(methodKey, ",M");
  auto structureOption = combineNamesForOptions(structureKey, ",x");
  auto gradientOption = combineNamesForOptions(calculateGradients, ",G");
  auto hessianOption = combineNamesForOptions(calculateHessian, ",H");
  auto atomicHessiansOption = combineNamesForOptions(calculateAtomicHessians, ",A");
  auto thermochemistryOption = combineNamesForOptions(thermochemistry, ",C");
  auto symmetryNumberOption = combineNamesForOptions(Utils::SettingsNames::symmetryNumber, ",S");
  auto temperatureOption = combineNamesForOptions(Utils::SettingsNames::temperature, ",T");
  auto suppressNormalModesOption = combineNamesForOptions(suppressNormalModesOutput, ",N");
  auto bondOrdersOption = combineNamesForOptions(calculateBondOrders, ",B");
  auto numberOrbitalSteersOption = combineNamesForOptions(numberOrbitalSteers, ",O");
  auto excitedStatesOption = combineNamesForOptions(calculateExcitedStates, ",E");
  auto numberOfExcitedStatesRootsOption = combineNamesForOptions(Utils::SettingsNames::numberOfEigenstates, ",r");
  auto initialGuessSpaceDimensionOption = combineNamesForOptions(Utils::SettingsNames::initialSubspaceDimension, ",n");
  auto excitedStatesSpinBlockOption = combineNamesForOptions(Utils::SettingsNames::spinBlock, ",b");
  auto excitedStatesDistanceThresholdOption = combineNamesForOptions(distanceThreshold, ",d");
  auto pruneBasisOption = combineNamesForOptions(Utils::SettingsNames::pruneBasis, ",y");
  auto maxMemoryOption = combineNamesForOptions(Utils::SettingsNames::maxMemory, ",g");
  auto energyThresholdOption = combineNamesForOptions(Utils::SettingsNames::energyThreshold, ",e");
  const auto* perturbativeThresholdOption = Utils::SettingsNames::perturbativeThreshold;
  auto wavefunctionOption = combineNamesForOptions(wavefunctionOutput, ",W");
  auto descriptionOption = combineNamesForOptions(calculationDescription, ",D");
  auto logOption = combineNamesForOptions(Utils::SettingsNames::loggerVerbosity, ",l");
  auto logFilenameOption = combineNamesForOptions(logFilename, ",f");

  // clang-format off
  pImpl_->desc_.add_options()
    (helpOption.c_str(), "prints this help message")
    (outputToFileOption.c_str(), "flag setting the saving of the output matrices as files")
    (chargeOption.c_str(), value<int>(), "sets the molecular charge")
    (multiplicityOption.c_str(), value<int>(), "sets the spin multiplicity")
    (unrestrictedOption.c_str(), "runs unrestricted calculation (UHF)")
    (parameterFileOption.c_str(), value<std::string>(), "sets the desired parameter path")
    (excitedStatesParameterFileOption.c_str(), value<std::string>(), "sets the desired parameter path")
    (mixerTypeOption.c_str(), value<std::string>(), "sets the desired SCF mixer for convergence acceleration")
    (maxScfIterationsOption.c_str(), value<int>(), "sets the maximum number of SCF iterations")
    (scfThresholdOption.c_str(), value<double>(), "sets the SCF convergence threshold for the energy difference (in hartree)")
    (Utils::SettingsNames::densityRmsdCriterion, value<double>(), "sets the SCF convergence threshold for the density matrix RMSD")
    (methodOption.c_str(), value<std::string>(), "sets the method used to perform the calculation")
    (structureOption.c_str(), value<std::string>()->value_name("structure"), "sets the path to the coordinates file")
    (gradientOption.c_str(), "activates the calculation of the atomic gradients")
    (hessianOption.c_str(), "activates the calculation of the molecular Hessian matrix")
    (atomicHessiansOption.c_str(), "activates the calculation of the atomic Hessian matrices")
    (thermochemistryOption.c_str(), "activates the calculation of thermochemical properties")
    (symmetryNumberOption.c_str(), value<int>(), "sets the molecular symmetry number for the thermochemical calculation")
    (temperatureOption.c_str(), value<double>(), "defines the temperature for the thermochemical calculation in Kelvin")
    (suppressNormalModesOption.c_str(), "suppresses the results of the hessian calculation as normal modes analysis")
    (bondOrdersOption.c_str(), "activates the calculation of the bond order matrix")
    (numberOrbitalSteersOption.c_str(), value<int>()->default_value(0), "sets the number of orbital steers to do, if 0, no steering")
    (numberOrbitalsToMix, value<int>(), "sets the number of orbital mixes to do per orbital steering")
    (numberOrbitalsToConsider, value<int>(), "sets the number of orbital around the fermi level that can be mixed")
    (maximalMixingAngle, value<double>(), "sets the maximal angle for the steering")
    (minimalMixingAngle, value<double>(), "sets the minimal angle for the steering")
    (maxIterationsSteering, value<int>(), "sets the maximal iteration number in the orbital steering run")
    (scfMixerSteering, value<std::string>(), "sets the SCF mixer (DIIS, eDIIS,...) in the orbital steering run")
    (excitedStatesOption.c_str(), "performs an excites state calculation with NDDO-CIS or TD-DFTB")
    (numberOfExcitedStatesRootsOption.c_str(), value<int>(), "sets the number of excited states to calculate")
    (initialGuessSpaceDimensionOption.c_str(), value<int>(), "sets the number of initial excited states guess vectors")
    (excitedStatesSpinBlockOption.c_str(), value<std::string>(), "sets excited states Hamiltonian spin-block to calculate (singlet/triplet/both)")
    (excitedStatesDistanceThresholdOption.c_str(), value<double>(), "sets the distance threshold after which the 2-e integrals contributions are neglected")
    (pruneBasisOption.c_str(), "sets whether the basis of singly excited determinants should be pruned")
    (energyThresholdOption.c_str(), value<double>(), "sets the threshold for pruning with an energy criterion in a.u.")
    (perturbativeThresholdOption, value<double>(), "in the pruning, sets the threshold for the perturbative energy contribution from secondary configurations")
    (maxMemoryOption.c_str(), value<double>(), "sets the maximum amount of memory that can be used by the calculation")
    (wavefunctionOption.c_str(), "outputs a molden file for the visualization of orbitals")
    (descriptionOption.c_str(), value<std::string>(), "sets a calculation description which will appear in the output")
    (logOption.c_str(), value<std::string>()->default_value("warning"), "sets whether warnings and errors are printed. Levels other than"
     "none, error, warning, output, debug throw an exception.")
    (logFilenameOption.c_str(), value<std::string>()->default_value(""), "Sets the name of the file where the logging is piped.");

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
    for (auto& c : result) {
      c = static_cast<char>(std::toupper(c));
    }
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
  if (pImpl_->vm_.count(Utils::SettingsNames::loggerVerbosity) > 0) {
    return pImpl_->vm_[Utils::SettingsNames::loggerVerbosity].as<std::string>();
  }
  return "";
}

std::string CommandLineOptions::getLogFilename() const {
  if (pImpl_->vm_.count(logFilename) > 0) {
    return pImpl_->vm_[Utils::SettingsNames::loggerVerbosity].as<std::string>();
  }
  return "";
}

int CommandLineOptions::getNumberOfOrbitalSteers() const {
  return pImpl_->vm_[numberOrbitalSteers].as<int>();
}

bool CommandLineOptions::gradientRequired() const {
  return pImpl_->vm_.count(calculateGradients) > 0;
}

bool CommandLineOptions::hessianRequired() const {
  return pImpl_->vm_.count(calculateHessian) > 0;
}

bool CommandLineOptions::atomicHessiansRequired() const {
  return pImpl_->vm_.count(calculateAtomicHessians) > 0;
}

bool CommandLineOptions::bondOrdersRequired() const {
  return pImpl_->vm_.count(calculateBondOrders) > 0;
}

bool CommandLineOptions::suppressNormalModes() const {
  return pImpl_->vm_.count(suppressNormalModesOutput) > 0;
}

bool CommandLineOptions::excitedStatesRequired() const {
  return pImpl_->vm_.count(calculateExcitedStates) > 0;
}

bool CommandLineOptions::orbitalSteeringRequired() const {
  return pImpl_->vm_[numberOrbitalSteers].as<int>() > 0;
}

bool CommandLineOptions::wavefunctionRequired() const {
  return pImpl_->vm_.count(wavefunctionOutput) > 0;
}

bool CommandLineOptions::thermochemistryRequired() const {
  return pImpl_->vm_.count(thermochemistry) > 0;
}

bool CommandLineOptions::pruneBasis() const {
  return pImpl_->vm_.count(Utils::SettingsNames::pruneBasis) > 0;
}

template<class CharPtrType, class StringType>
std::string CommandLineOptions::combineNamesForOptions(CharPtrType nameOfOption, StringType abbreviatedOption) const {
  auto combinedString = static_cast<std::string>(nameOfOption) + static_cast<std::string>(abbreviatedOption);
  return combinedString;
}

void CommandLineOptions::updateLogger(Core::Log& log) const {
  auto verbosity = getLoggerVerbosity();
  auto filename = getLogFilename();
  enum Type { Error, Output };

  auto getSink = [&](Type type) -> Core::Log::SinkPtr {
    if (filename.empty()) {
      return type == Type::Error ? Core::Log::cerrSink() : Core::Log::coutSink();
    }
    return Core::Log::fileSink(filename);
  };

  if (verbosity.empty()) {
    return;
  }

  log = Core::Log::silent();
  if (verbosity == "none") {
    return;
  }
  if (verbosity == "error") {
    log.error.add("cerr", getSink(Type::Error));
  }
  else if (verbosity == "warning") {
    log.error.add("cerr", getSink(Type::Error));
    log.warning.add("cout", getSink(Type::Output));
  }
  else if (verbosity == "output") {
    log.error.add("cerr", getSink(Type::Error));
    log.warning.add("cout", getSink(Type::Output));
    log.output.add("cout", getSink(Type::Output));
  }
  else if (verbosity == "debug") {
    log.error.add("cerr", getSink(Type::Error));
    log.warning.add("cout", getSink(Type::Output));
    log.output.add("cout", getSink(Type::Output));
    log.debug.add("cout", getSink(Type::Output));
  }
  else {
    throw std::runtime_error("Logging level '" + verbosity + "' does not exist.");
  }
}

void CommandLineOptions::updateSettings(Utils::Settings& settingsToUpdate) const {
  using namespace Utils::SettingsNames;
  if (validOptionToSet(molecularCharge, settingsToUpdate)) {
    settingsToUpdate.modifyInt(molecularCharge, pImpl_->vm_[molecularCharge].as<int>());
  }
  if (validOptionToSet(spinMultiplicity, settingsToUpdate)) {
    settingsToUpdate.modifyInt(spinMultiplicity, pImpl_->vm_[spinMultiplicity].as<int>());
  }
  if (settingsToUpdate.valueExists(spinMode)) {
    if (pImpl_->vm_.count(unrestrictedCalculation) > 0) {
      settingsToUpdate.modifyString(spinMode, Utils::SpinModeInterpreter::getStringFromSpinMode(Utils::SpinMode::Unrestricted));
    }
    else {
      settingsToUpdate.modifyString(spinMode, Utils::SpinModeInterpreter::getStringFromSpinMode(Utils::SpinMode::Restricted));
    }
  }
  if (validOptionToSet(mixer, settingsToUpdate)) {
    settingsToUpdate.modifyString(mixer, pImpl_->vm_[mixer].as<std::string>());
  }
  if (validOptionToSet(maxScfIterations, settingsToUpdate)) {
    settingsToUpdate.modifyInt(maxScfIterations, pImpl_->vm_[maxScfIterations].as<int>());
  }
  if (validOptionToSet(selfConsistenceCriterion, settingsToUpdate)) {
    settingsToUpdate.modifyDouble(selfConsistenceCriterion, pImpl_->vm_[selfConsistenceCriterion].as<double>());
  }
  if (validOptionToSet(densityRmsdCriterion, settingsToUpdate)) {
    settingsToUpdate.modifyDouble(densityRmsdCriterion, pImpl_->vm_[densityRmsdCriterion].as<double>());
  }
  if (validOptionToSet(methodParameters, settingsToUpdate)) {
    settingsToUpdate.modifyString(methodParameters, pImpl_->vm_[methodParameters].as<std::string>());
  }
  if (validOptionToSet(temperature, settingsToUpdate)) {
    settingsToUpdate.modifyDouble(temperature, pImpl_->vm_[temperature].as<double>());
  }
  if (validOptionToSet(symmetryNumber, settingsToUpdate)) {
    settingsToUpdate.modifyInt(symmetryNumber, pImpl_->vm_[symmetryNumber].as<int>());
  }
}

void CommandLineOptions::updateExcitedStatesSettings(Utils::Settings& settingsToUpdate) const {
  using namespace Utils;
  if (validOptionToSet(SettingsNames::excitedStatesParamFile, settingsToUpdate)) {
    settingsToUpdate.modifyString(SettingsNames::excitedStatesParamFile,
                                  pImpl_->vm_[SettingsNames::excitedStatesParamFile].as<std::string>());
  }
  if (validOptionToSet(SettingsNames::numberOfEigenstates, settingsToUpdate)) {
    settingsToUpdate.modifyInt(SettingsNames::numberOfEigenstates, pImpl_->vm_[SettingsNames::numberOfEigenstates].as<int>());
  }
  if (validOptionToSet(SettingsNames::initialSubspaceDimension, settingsToUpdate)) {
    settingsToUpdate.modifyInt(SettingsNames::initialSubspaceDimension,
                               pImpl_->vm_[SettingsNames::initialSubspaceDimension].as<int>());
  }
  if (validOptionToSet(SettingsNames::spinBlock, settingsToUpdate)) {
    settingsToUpdate.modifyString(SettingsNames::spinBlock, pImpl_->vm_[SettingsNames::spinBlock].as<std::string>());
  }
  if (validOptionToSet(SettingsNames::energyThreshold, settingsToUpdate)) {
    settingsToUpdate.modifyDouble(SettingsNames::energyThreshold, pImpl_->vm_[SettingsNames::energyThreshold].as<double>());
  }
  if (validOptionToSet(SettingsNames::perturbativeThreshold, settingsToUpdate)) {
    settingsToUpdate.modifyDouble(SettingsNames::perturbativeThreshold,
                                  pImpl_->vm_[SettingsNames::perturbativeThreshold].as<double>());
  }
  if (validOptionToSet(SettingsNames::maxMemory, settingsToUpdate)) {
    settingsToUpdate.modifyDouble(SettingsNames::maxMemory, pImpl_->vm_[SettingsNames::maxMemory].as<double>());
  }
  if (validOptionToSet(distanceThreshold, settingsToUpdate)) {
    settingsToUpdate.modifyDouble(distanceThreshold, pImpl_->vm_[distanceThreshold].as<double>());
  }
  if (pruneBasis()) { // If new pruning methods are implemented, this needs to change.
    if (validOptionToSet(SettingsNames::pruneBasis, settingsToUpdate)) {
      settingsToUpdate.modifyString(SettingsNames::pruneBasis, "energy");
    }
    else {
      throw std::runtime_error("Basis pruning only available with TD-DFTB methods");
    }
  }
}

void CommandLineOptions::updateOrbitalSteeringSettings(Utils::Settings& settingsToUpdate) const {
  using namespace Utils;
  if (validOptionToSet(maximalMixingAngle, settingsToUpdate)) {
    settingsToUpdate.modifyDouble(maximalMixingAngle, pImpl_->vm_[maximalMixingAngle].as<double>());
  }
  if (validOptionToSet(minimalMixingAngle, settingsToUpdate)) {
    settingsToUpdate.modifyDouble(minimalMixingAngle, pImpl_->vm_[minimalMixingAngle].as<double>());
  }
  if (validOptionToSet(numberOrbitalsToMix, settingsToUpdate)) {
    settingsToUpdate.modifyInt(numberOrbitalsToMix, pImpl_->vm_[numberOrbitalsToMix].as<int>());
  }
  if (validOptionToSet(numberOrbitalsToConsider, settingsToUpdate)) {
    settingsToUpdate.modifyInt(numberOrbitalsToConsider, pImpl_->vm_[numberOrbitalsToConsider].as<int>());
  }
  if (validOptionToSet(maxIterationsSteering, settingsToUpdate)) {
    settingsToUpdate.modifyInt(maxIterationsSteering, pImpl_->vm_[maxIterationsSteering].as<int>());
  }
  if (validOptionToSet(scfMixerSteering, settingsToUpdate)) {
    using namespace Utils::SettingsNames;
    std::vector<std::string> availableMixers = {ScfMixers::noMixer, ScfMixers::diis, ScfMixers::ediis, ScfMixers::ediisDiis};
    std::string caseInsensitiveMixer = scfMixerSteering;
    for (auto& mixer : availableMixers) {
      std::transform(mixer.begin(), mixer.end(), mixer.begin(), ::toupper);
    }
    std::transform(caseInsensitiveMixer.begin(), caseInsensitiveMixer.end(), caseInsensitiveMixer.begin(), ::toupper);
    auto it = std::find(availableMixers.begin(), availableMixers.end(), caseInsensitiveMixer);
    if (it == availableMixers.end()) {
      throw std::runtime_error("No mixer called " + std::string(scfMixerSteering) + " available.");
    }
    settingsToUpdate.modifyString(scfMixerSteering, pImpl_->vm_[scfMixerSteering].as<std::string>());
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
  std::string callStatement;
  for (int i = 0; i < argc; ++i) {
    callStatement += " " + std::string(argv[i]);
  }
  return callStatement;
}

} // namespace Sparrow
} // namespace Scine
