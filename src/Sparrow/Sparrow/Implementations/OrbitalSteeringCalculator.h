/**
 * @file Calculator.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_ORBITALSTEERINGCALCULATOR_H
#define SPARROW_ORBITALSTEERINGCALCULATOR_H

#include <Core/Interfaces/CalculatorWithReference.h>

namespace Scine {
namespace Utils {
class Results;
} // namespace Utils
namespace Sparrow {

class GenericMethodWrapper;
class OrbitalSteeringSettings;

class OrbitalSteeringCalculator final : public Core::CalculatorWithReference {
 public:
  struct InvalidCalculatorTypeForOrbitalSteerer final : public std::exception {
    InvalidCalculatorTypeForOrbitalSteerer();
    const char* what() const noexcept final {
      return error.c_str();
    }
    std::string error;
  };
  static constexpr const char* model = "orbital_steering";
  /// @brief Default constructor.
  OrbitalSteeringCalculator();
  /// @brief Virtual destructor.
  ~OrbitalSteeringCalculator() final;
  /**
   * @brief Sets the calculator to be used to perform the reference calculation.
   * In the derived classes care must be taken that the case where a method does not accept
   * some calculator types (i.e CIS with DFT, or TDDFT with HF) is checked and handled.
   * @throws if the referenceCalculator is not a valid reference calculator for this class instance.
   */
  void setReferenceCalculator(std::shared_ptr<Core::Calculator> referenceCalculator) final;
  /**
   * @brief Performs a reference calculation.
   */
  void referenceCalculation() final;

  /**
   * @brief Accessor for the reference calculator.
   * @return Core::Calculator& The reference calculator.
   */
  Core::Calculator& getReferenceCalculator() final;
  /**
   * @brief Constant accessor for the reference calculator.
   * @return const Core::Calculator& The reference calculator.
   */
  const Core::Calculator& getReferenceCalculator() const final;

  /**
   * @brief The main function running the calculation with reference.
   * @returns A const-ref of stored (and newly calculated) Results.
   */
  const Utils::Results& calculate() final;

  /**
   * @brief Getter for the name of the calculator with reference.
   * @return Returns the name of the calculator with reference.
   */
  std::string name() const final;

  /**
   * @brief Accessor for the settings.
   * @return Utils::Settings& The settings.
   */
  Utils::Settings& settings() final;
  /**
   * @brief Constant accessor for the settings.
   * @return const Utils::Settings& The settings.
   */
  const Utils::Settings& settings() const final;
  /**
   * @brief Method to apply the settings stored in the settings data structure.
   */
  void applySettings() final;
  /**
   * @brief Accessor for the saved instance of Utils::Results.
   * @return Utils::Results& The results of the previous calculation.
   */
  Utils::Results& results() final;
  /**
   * @brief Constant accessor for the Utils::Results.
   * @return const Utils::Results& The results of the previous calculation.
   */
  const Utils::Results& results() const final;

 private:
  void logSteering(double oldEnergy, double newEnergy);
  std::shared_ptr<GenericMethodWrapper> method_;
  std::unique_ptr<OrbitalSteeringSettings> settings_;
  int numberOfCalculations_;
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_ORBITALSTEERINGCALCULATOR_H
