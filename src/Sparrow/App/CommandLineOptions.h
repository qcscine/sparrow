/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_COMMANDLINEOPTIONS_H
#define SPARROW_COMMANDLINEOPTIONS_H

#include <memory>
#include <ostream>

namespace Scine {
namespace Utils {
class Settings;
} // namespace Utils
namespace Sparrow {

/**
 * @brief Class to parse the command line options for non-default options and passes them to a Util::Settings class.
 * This class uses the pImpl idiom to hide the boost::program_options dependency.
 */
class CommandLineOptions {
 public:
  /**
   * @brief Class constructor, parses the command line arguments and maps them to the according setting.
   * @param argv the vector with the argument strings.
   * @param argc the number of arguments.
   */
  CommandLineOptions(int argc, char* argv[]);
  ~CommandLineOptions();

  /** @brief returns the command call used to run the program. */
  std::string getCallStatement() const;
  /** @brief returns whether the help flag option has been set. */
  bool helpRequired() const;
  /** @brief returns whether the matrices should be saved as files. */
  bool outputToFileRequired() const;

  /** @brief returns the method name given as command-line argument. */
  std::string getSelectedMethodName() const;
  /** @brief returns the xyz file containing the coordinates with the desired structures. */
  std::string getStructureCoordinatesFile() const;
  /** @brief returns the desired calculation description. */
  std::string getCalculationDescription() const;
  /** @brief returns the desired verbosity of the logging. */
  std::string getLoggerVerbosity() const;
  /** @brief returns whether the gradients have to be computed. */
  bool gradientRequired() const;
  /** @brief returns whether the hessian matrix has to be computed. */
  bool hessianRequired() const;
  /** @brief returns whether the bond order matrix has to be computed. */
  bool bondOrdersRequired() const;
  /** @brief returns whether the normal modes output is printed or suppressed. */
  bool suppressNormalModes() const;
  /** @brief returns whether the wavefunction output is printed as a molden file. */
  bool wavefunctionRequired() const;
  /** @brief returns whether the thermochemical properties are calculated. */
  bool thermochemistryRequired() const;

  /** @brief updates a setting with the option parsed from the command line. */
  void updateSettings(Utils::Settings& settingsToUpdate) const;

  /** @brief prints the help message. */
  void printHelp(std::ostream& out) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> pImpl_;

  /// @brief Parses the command line to generate a call statement.
  std::string generateCallStatement(int argc, char* argv[]) const;

  /// templated function to allow using with string and const char pointers as argument.
  /// Checks whether the option given by optionIdentifier can be used in the corresponding settings.
  template<class CharType>
  bool validOptionToSet(CharType optionIdentifier, const Utils::Settings& settings) const;

  /// Combines two character types to form a single const char*.
  template<class CharPtrType, class StringType>
  std::string combineNamesForOptions(CharPtrType nameOfOption, StringType abbreviatedOption) const;
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_COMMANDLINEOPTIONS_H
