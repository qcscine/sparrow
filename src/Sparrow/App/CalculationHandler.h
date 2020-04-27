/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_CALCULATIONHANDLER_H
#define SPARROW_CALCULATIONHANDLER_H

#include <Utils/CalculatorBasics/Results.h>
#include <chrono>
#include <exception>
#include <memory>
#include <ostream>

namespace Scine {
namespace Utils {
class Settings;
} // namespace Utils

namespace Core {
class Calculator;
} // namespace Core
namespace Sparrow {

class SparrowInitializer;
class CommandLineOptions;

/**
 * @brief Exception thrown if a non available method is requested.
 */
class MethodNotAvailableException : public std::exception {
 public:
  explicit MethodNotAvailableException(std::string methodName)
    : methodName_("Method " + std::move(methodName) + " is not available.\n" + "\nSparrow flies away...") {
  }

  const char* what() const noexcept final {
    return methodName_.c_str();
  }

 private:
  std::string methodName_;
};

struct FileInaccessibleException final : public std::exception {
 public:
  explicit FileInaccessibleException(std::string fileName)
    : description_("File " + std::move(fileName) + " non accessible.") {
  }
  const char* what() const noexcept final {
    return description_.c_str();
  }

 private:
  std::string description_;
};
/*
 * @class CalculationHandler @file CalculationHandler.h
 * @brief Class handling the main calculation routines.
 */
class CalculationHandler {
 public:
  //! @brief deleted default constructor;
  CalculationHandler() = delete;
  /**
   * @brief This constructor initializes the desired method with the settings
   *         read from the command line.
   */
  CalculationHandler(CommandLineOptions& options, SparrowInitializer& initializer);
  /**
   * @brief start a calculation with the options and description read by the command line.
   */
  void calculate(std::ostream& out);

 private:
  void assignPropertiesToCalculate();
  void assignSettings();
  void printSettings(std::ostream& out, const Utils::Settings& settings, std::string commentChar = "") const;
  void printHeader(std::ostream& out, std::string commenChar = "") const;
  void printFooter(std::ostream& out) const;
  void printTime(std::ostream& out) const;
  void printResultsToFile() const;
  void printPrettyResults(std::ostream& out) const;
  void printFrequencyAnalysis(std::ostream& out, const Utils::HessianMatrix& matrix) const;
  void printWavefunction() const;
  CommandLineOptions& commandLineOptions_;
  std::shared_ptr<Core::Calculator> methodWrapper_;
  Utils::Results results_;
  std::chrono::time_point<std::chrono::system_clock> start_, end_;
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_CALCULATIONHANDLER_H
