/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_CISLINEARRESPONSETIMEDEPENDENTCALCULATOR_H
#define SPARROW_CISLINEARRESPONSETIMEDEPENDENTCALCULATOR_H

#include "CISSigmaVectorEvaluator.h"
#include <Core/Interfaces/CalculatorWithReference.h>
#include <Sparrow/Implementations/TimeDependent/LinearResponseCalculator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Math/IterativeDiagonalizer/SpinAdaptedEigenContainer.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <chrono>
#include <memory>
#include <numeric>
namespace Scine {
namespace Core {
class Calculator;
} // namespace Core
namespace Utils {
class Settings;
class Results;
namespace LcaoUtils {
class ElectronicOccupation;
}
enum class SpinTransition;
} // namespace Utils
namespace Sparrow {

class NDDOMethodWrapper;

/**
 * @class CISLinearResponseTimeDependentCalculator
 * @file CISLinearResponseTimeDependentCalculator.h
 * @brief Top class responsible for CIS calculations of excited states.
 * This class calculates the CIS excited states via a direct algorithm to avoid
 * the O(N^5) 4-index transformation.
 * This class works for all the NDDO methods implemented in Sparrow: MNDO, AM1, PM3, RM1, PM6.
 * @pre The Sparrow::NDDOMethodWrapper given in the constructor must contain the calculated 2-electron integrals.
 * @pre The Core::Calculator given as reference must derive from Sparrow::NDDOMethodWrapper,
 *      so it still does not work with any Core::Calculator.
 */
class CISLinearResponseTimeDependentCalculator final : public LinearResponseCalculator {
 public:
  constexpr static const char* model = "CIS-NDDO";
  /// TODO: Active space?
  CISLinearResponseTimeDependentCalculator();
  ~CISLinearResponseTimeDependentCalculator() final = default;

  /**
   * @brief Accessor for the ground state calculator.
   * @return Core::Calculator& The reference calculator.
   */
  Core::Calculator& getReferenceCalculator() final;
  /**
   * @brief Constant accessor for the ground state calculator.
   * @return const Core::Calculator& The reference calculator.
   */
  const Core::Calculator& getReferenceCalculator() const final;
  /**
   * @brief Sets the ground state method to calculate excited states.
   * @param method A NDDO method.
   * If the method is not NDDO, an instance InvalidCalculatorTypeForRHFCIS is thrown.
   */
  void setReferenceCalculator(std::shared_ptr<Core::Calculator> method) final;
  /**
   * @brief This function gives the chance to perform a reference calculation.
   * This is not needed if the calculator has already performed a RHF calculation.
   * @pre nddoMethod_ must already be initialized and equipped with a structure to calculate, i.e.
   *      the function setStructure(Utils::AtomCollection) must already be called.
   */
  void referenceCalculation() final;

  /**
   * @brief Accessor for the Results stored in this method wrapper.
   * @returns Utils::Results& The results of the previous calculation.
   */
  Utils::Results& results() final;
  /**
   * @brief Const accessor for the Results used in this method wrapper.
   * @returns const Utils::Results& The results of the previous calculation.
   */
  const Utils::Results& results() const final;
  /**
   * @brief Accessor for the underlying settings.
   */
  Utils::Settings& settings() final;
  /**
   * @brief Constant accessor for the underlying settings.
   */
  const Utils::Settings& settings() const final;
  /**
   * @brief Getter for the model name, in this case "CIS".
   */
  std::string name() const final;

  /**
   * @brief Apply the settings contained in the settings_ member.
   */
  void applySettings() final;

  /**
   * @brief Solves the first roots of the CIS Matrix in a direct way.
   * @param numberOfEnergyLevels The desired amount of energy levels to compute. If equal to 0
   *        then all the energy levels are calculated.
   * @return A pair of a Eigen::VectorXd containing the eigenvalues and a
   *         correspondingly sorted Eigen::MatrixXd containing the eigenvectors.
   * This function calculates the first roots of the CIS matrix in a direct way,
   * i.e. directly from the AO integrals without any storage of the CIS matrix or
   * 2-index transformation of the AO integrals to MO integrals.
   */
  template<Utils::Reference restrictedness>
  Utils::ElectronicTransitionResult solve(Utils::SpinAdaptedContainer<restrictedness, Eigen::VectorXd> energyDifferenceVector,
                                          int numberOfEnergyLevels = 0, int initialSubspaceDimension = 0,
                                          Utils::SpinTransition spinBlock = Utils::SpinTransition::Singlet);

  /**
   * @brief Overridden calculate method, inherited from Core::CalculatorWithReference.
   * This method works through the settings. It is equivalent to a non-polymorphically
   * available solve(int, int, Utils::SpinTransition) call.
   * @return Returns a Utils::Results class containing a Utils::SpinAdaptedEigenContainer
   *         instance.
   */
  const Utils::Results& calculate() final;

  /**
   * @brief Sets the guess to be used in the next calculation. If empty, diagonal dominant guess will be used.
   */
  void setGuess(std::shared_ptr<GuessSpecifier> guessVectorMatrix) final;
  /**
   * @brief Returns the guess in the full singles space (no pruning)
   */
  auto getGuess() const -> std::shared_ptr<GuessSpecifier> final;

 private:
  void setExcitedStatesParam(Utils::Reference restrictedness, Utils::SpinTransition spinBlock);
  void prepareIntegralScreening();
  void checkMemoryRequirement(int excitationsDim, int numberOfEnergyLevels);
  template<Utils::Reference restrictedness>
  void generateTransitionDipoleMoments(Utils::ElectronicTransitionResult& excitedStatesResults, const CISData& cisData,
                                       Utils::SpinTransition spinBlock) const;

  std::shared_ptr<NDDOMethodWrapper> nddoMethod_;
  std::unique_ptr<Utils::Settings> settings_;
  std::unique_ptr<CISData> cisData_;
  std::shared_ptr<GuessSpecifier> guess_;
  std::vector<std::multimap<double, int, std::greater<double>>> integralsThresholds_;
  std::vector<int> orderMap_;
  ExcitedStatesParam excitedStatesParam_{1., 1., 1.};
  Utils::Results results_;
};
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_CISLINEARRESPONSETIMEDEPENDENTCALCULATOR_H
