/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_TDDFTBCALCULATOR_H
#define SPARROW_TDDFTBCALCULATOR_H

#include <Core/Interfaces/CalculatorWithReference.h>
#include <Sparrow/Implementations/TimeDependent/LinearResponseCalculator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Math/IterativeDiagonalizer/SigmaVectorEvaluator.h>
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
class Excitation;
class DipoleMatrix;
namespace LcaoUtils {
class ElectronicOccupation;
}
enum class SpinTransition;
} // namespace Utils
namespace Sparrow {

class DFTBMethodWrapper;
class TDDFTBData;
template<Utils::Reference restrictedness>
class OrderedInput;
template<Utils::Reference restrictedness>
class TDDFTBEigenvalueSolver;

/**
 * @class TDDFTBCalculator
 * @file TDDFTBCalculator.h
 * @brief Top class responsible for TDDFTB calculations of excited states.
 *
 * Implementation details for the singlet closed-shell case:
 * R. Rüger, E. van Lenthe, Y. Lu, J. Frenzel, T. Heine, L. Visscher,
 * Efficient Calculation of Electronic Absorption Spectra by Means of Intensity-Selected
 * Time-Dependent Density Functional Tight Binding, JCTC, 2015,
 * https://doi.org/10.1021/ct500838h
 * Adaptations to Triplet and Unrestricted cases:
 * T. A. Niehaus, S. Suhai, F. Della Sala, P. Lugli, M. Elstner, G. Seifert, and Th. Frauenheim,
 * Tight-binding approach to time-dependent density-functional response theory,
 * Phys. Rev. B 63, 085108, 2001,
 * https://doi.org/10.1103/PhysRevB.63.085108
 *
 * This class works for all the DFTB methods implemented in Sparrow: DFTB0, DFTB2, DFTB3.
 *
 * Excited states in DFTB0 are just the energy difference between the orbitals.
 * The DFTB3-specific corrections are neglected, so that the excited states are calculated just at
 * te SCC-DFTB2 level. This is justified by
 * Y. Nishimoto, Time-dependent density-functional tight-binding method
 * with the third-order expansion of electron density, J. Chem. Phys. 143, 094108 (2015).
 *
 * @pre The Sparrow::DFTBMethodWrapper given in the constructor.
 * @pre The Core::Calculator given as reference must derive from Sparrow::DFTBMethodWrapper.
 *      Failure results in a throw.
 */
class TDDFTBCalculator final : public LinearResponseCalculator {
 public:
  constexpr static const char* model = "TD-DFTB";
  TDDFTBCalculator();
  ~TDDFTBCalculator() final;

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
   * @param method A DFTB method.
   * If the method is not DFTB, an instance InvalidCalculatorTypeForTDDFTB is thrown.
   */
  void setReferenceCalculator(std::shared_ptr<Core::Calculator> method) final;
  /**
   * @brief This function gives the chance to perform a reference calculation.
   * This is not needed if the calculator has already performed a RHF calculation.
   * @pre dftbMethod_ must already be initialized and equipped with a structure to calculate, i.e.
   *      the function setStructure(Utils::AtomCollection) must already be called.
   */
  void referenceCalculation() final;

  /**
   * @brief Accessor for the Results stored.
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
   * @brief Getter for the model name, in this case "TDDFTB".
   */
  std::string name() const final;

  /**
   * @brief Apply the settings contained in the settings_ member.
   */
  void applySettings() final;

  template<Utils::Reference restrictedness>
  Utils::ElectronicTransitionResult solve(Utils::SpinAdaptedContainer<restrictedness, Eigen::VectorXd> energyDifferenceVector,
                                          int numberOfEnergyLevels = 0, int initialSubspaceDimension = 0,
                                          Utils::SpinTransition spinBlock = Utils::SpinTransition::Singlet);

  /**
   * @brief Overridden calculate method, inherited from Core::CalculatorWithReference.
   * This method works through the settings. It is equivalent to a non-polymorphically
   * available
   * solve(Utils::SpinAdaptedContainer<restrictedness, Eigen::VectorXd>, int, int, Utils::SpinTransition)
   * call.
   * @return Returns a Utils::Results class containing a Utils::SpinAdaptedEigenContainer
   *         instance.
   */
  const Utils::Results& calculate() final;

  /**
   * @brief Sets the guess to be used in the next calculation. If empty, diagonal dominant guess will be used.
   */
  void setGuess(std::shared_ptr<GuessSpecifier> guessVectorMatrix) final;
  auto getGuess() const -> std::shared_ptr<GuessSpecifier> final;

 private:
  void checkMemoryRequirement(int excitationsDim, int numberOfEnergyLevels);
  template<Utils::Reference restrictedness>
  void generateTransitionDipoleMoments(Utils::ElectronicTransitionResult& excitedStatesResults,
                                       const std::vector<Utils::Excitation>& excitations) const;
  auto isDFTB0(std::shared_ptr<DFTBMethodWrapper> method) const -> bool;

  template<Utils::Reference restrictedness>
  void evalImpl(Utils::SpinAdaptedElectronicTransitionResult& transitionResult,
                const Eigen::MatrixXd& transitionCharges, std::shared_ptr<Utils::DipoleMatrix> dipoleMatrix);
  template<Utils::Reference restrictedness>
  void solveEigenvalueProblem(Utils::SpinAdaptedElectronicTransitionResult& transitionResult,
                              const std::unique_ptr<TDDFTBEigenvalueSolver<restrictedness>>& solver, int numberOfRoots,
                              int initialSubspaceDimension);
  std::shared_ptr<DFTBMethodWrapper> dftbMethod_;
  std::shared_ptr<GuessSpecifier> guess_;
  std::unique_ptr<Utils::Settings> settings_;
  std::unique_ptr<TDDFTBData> tddftbData_;
  std::vector<int> orderMap_;
  std::vector<Utils::Excitation> excitations_;
  Utils::Results results_;
};
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_TDDFTBCALCULATOR_H
