/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_TDDFTBEIGENVALUESOLVER_H
#define SPARROW_TDDFTBEIGENVALUESOLVER_H

#include "BasisPruner.h"
#include "TDDFTBSigmaVectorEvaluator.h"
#include <Sparrow/Implementations/TimeDependent/DiagonalPreconditionerEvaluator.h>
#include <Sparrow/Implementations/TimeDependent/LinearResponseCalculator.h>
#include <Sparrow/Implementations/TimeDependent/LinearResponseSettings.h>
#include <Sparrow/Implementations/TimeDependent/TimeDependentUtils.h>
#include <Utils/Math/IterativeDiagonalizer/DavidsonDiagonalizer.h>
#include <Utils/Math/IterativeDiagonalizer/DiagonalizerSettings.h>

namespace Scine {
namespace Sparrow {

template<Utils::Reference restrictedness>
class TDDFTBEigenvalueSolver {
 public:
  TDDFTBEigenvalueSolver(const Utils::Settings& settings, std::shared_ptr<Eigen::MatrixXd> gammaMatrix,
                         std::shared_ptr<Eigen::VectorXd> spinConstants, const OrderedInput<restrictedness>& orderedInput,
                         std::shared_ptr<Utils::DipoleMatrix> dipoleMatrixMO, Core::Log& log)
    : settings_(settings),
      gammaMatrix_(std::move(gammaMatrix)),
      spinConstants_(std::move(spinConstants)),
      dipoleMatrix_(std::move(dipoleMatrixMO)),
      input_(orderedInput),
      log_(log) {
  }

  /**
   * @brief Solves the first roots of the TDDFTB Matrix without any pruning.
   * @param numberOfEnergyLevels The desired amount of energy levels to compute. If equal to 0
   *        then all the energy levels are calculated.
   * @return A pair of a Eigen::VectorXd containing the eigenvalues and a
   *         correspondingly sorted Eigen::MatrixXd containing the eigenvectors.
   * This function calculates the first roots of the TDDFTB matrix
   */
  auto solve(int numberOfEnergyLevels, int initialSubspaceDimension,
             Utils::SpinTransition spinBlock = Utils::SpinTransition::Singlet) -> Utils::ElectronicTransitionResult;

  auto solveDftb0(int numberOfRoots) const -> Utils::ElectronicTransitionResult;

  void setGuess(std::shared_ptr<LinearResponseCalculator::GuessSpecifier> guess) {
    guess_ = std::move(guess);
  }

 private:
  void checkAndCorrectNumberOfRoots(int& numberOfEnergyLevels, int& initialSubspaceDimension, int nConfigurations);
  void generateTransitionDipoleMoments(Utils::ElectronicTransitionResult& excitedStatesResults,
                                       Utils::SpinTransition spinBlock) const;
  /**
   * @brief Transforms the results of the diagonalization of the \Omega matrix to the electronic transition results
   * In particular:
   * - Transition energies:
   * w_I = sqrt(v_I)
   * w_I: I-th electronic vertical transition energy
   * v_i: I-th eigenvalue
   *
   * - Coefficients:
   * c_{ia}^I = sqrt((e_a - e_i) / w_I) * F_{ia}^I
   * c_{ia}^I: coefficient of the \Phi_ia determinant in the I-th excited state
   * e_a - e_i: single particle energy difference for the \Phi_ia determinant
   * w_I: I-th electronic vertical transition energy
   * F_{ia}^I: ia-th coefficient of the I-th eigenvector
   * The c_{ia}^I coefficients are also used to compute the transition dipole moment.
   */
  Utils::ElectronicTransitionResult formExcitedStatesResults(const Utils::EigenContainer& resultsFromDavidson,
                                                             Utils::SpinTransition spinBlock) const;

  auto generateGuess(int nConfigurations, Utils::SpinTransition spinBlock = Utils::SpinTransition::Singlet)
      -> boost::optional<Eigen::MatrixXd>;
  const Utils::Settings& settings_;
  std::shared_ptr<Eigen::MatrixXd> gammaMatrix_;
  std::shared_ptr<Eigen::VectorXd> spinConstants_;
  std::shared_ptr<Utils::DipoleMatrix> dipoleMatrix_;
  std::shared_ptr<LinearResponseCalculator::GuessSpecifier> guess_;
  OrderedInput<restrictedness> input_;
  Core::Log& log_;
};

template<Utils::Reference restrictedness>
Utils::ElectronicTransitionResult TDDFTBEigenvalueSolver<restrictedness>::solve(int numberOfEnergyLevels,
                                                                                int initialSubspaceDimension,
                                                                                Utils::SpinTransition spinBlock) {
  const int nConfigurations = input_.energyDifferences().size();
  TDDFTBType type = settings_.getBool("tda") ? TDDFTBType::TDA : TDDFTBType::TDDFTB;

  checkAndCorrectNumberOfRoots(numberOfEnergyLevels, initialSubspaceDimension, nConfigurations);

  Utils::NonOrthogonalDavidson diagonalizer(numberOfEnergyLevels, nConfigurations);
  diagonalizer.settings().modifyInt(Utils::initialGuessDimensionOption, initialSubspaceDimension);
  diagonalizer.settings().modifyDouble(Utils::residualNormToleranceOption, settings_.getDouble(convergence));
  diagonalizer.settings().modifyString("gep_algo", settings_.getString("gep_algo"));
  if (settings_.getInt(Utils::SettingsNames::maxDavidsonIterations) != 0) {
    diagonalizer.settings().modifyInt(Utils::SettingsNames::maxDavidsonIterations,
                                      settings_.getInt(Utils::SettingsNames::maxDavidsonIterations));
  }

  diagonalizer.setGuess(generateGuess(nConfigurations, spinBlock));

  diagonalizer.setPreconditionerEvaluator(std::make_unique<DiagonalPreconditionerEvaluator>(input_.energyDifferences()));

  diagonalizer.setSigmaVectorEvaluator(
      std::make_unique<TDDFTBSigmaVectorEvaluator<restrictedness>>(gammaMatrix_, spinConstants_, input_, spinBlock, type));

  auto eigenvalueProblemResult = diagonalizer.solve(log_);

  auto excitedStates = formExcitedStatesResults(eigenvalueProblemResult, spinBlock);

  return excitedStates;
}

template<Utils::Reference restrictedness>
Utils::ElectronicTransitionResult TDDFTBEigenvalueSolver<restrictedness>::solveDftb0(int numberOfEnergyLevels) const {
  const int nRoots = std::min(numberOfEnergyLevels, int(input_.energyDifferences().size()));
  Utils::ElectronicTransitionResult result;
  Utils::EigenContainer diagonalResult;
  result.eigenStates.eigenValues = input_.energyDifferences().head(nRoots);
  result.eigenStates.eigenVectors = Eigen::MatrixXd::Identity(input_.energyDifferences().size(), nRoots);

  // DFTB0 has no spin-unrestricted, and eigenvalues correspond to the identity matrix (no coupling)
  // so singlet and triplets have the same energy. Oscillator strength printed for singlet. 0 for triplet.
  generateTransitionDipoleMoments(result, Utils::SpinTransition::Singlet);

  return result;
}

template<Utils::Reference restrictedness>
void TDDFTBEigenvalueSolver<restrictedness>::checkAndCorrectNumberOfRoots(int& numberOfEnergyLevels,
                                                                          int& initialSubspaceDimension, int nConfigurations) {
  // If 0 is given, then give all the energy levels
  if (numberOfEnergyLevels == 0 || numberOfEnergyLevels > nConfigurations) {
    numberOfEnergyLevels = nConfigurations;
  }
  if (initialSubspaceDimension == 0 || initialSubspaceDimension < numberOfEnergyLevels) {
    initialSubspaceDimension = numberOfEnergyLevels;
  }
  if (initialSubspaceDimension > nConfigurations) {
    initialSubspaceDimension = nConfigurations;
  }
}

template<Utils::Reference restrictedness>
Utils::ElectronicTransitionResult
TDDFTBEigenvalueSolver<restrictedness>::formExcitedStatesResults(const Utils::EigenContainer& resultsFromDavidson,
                                                                 Utils::SpinTransition spinBlock) const {
  Utils::ElectronicTransitionResult excitedStates;

  if (settings_.getBool("tda")) {
    excitedStates.eigenStates.eigenVectors = resultsFromDavidson.eigenVectors;
    excitedStates.eigenStates.eigenValues = resultsFromDavidson.eigenValues;
  }
  else {
    // c_{ia}^I = sqrt((e_a - e_i) / w_I) * F_{ia}^I
    excitedStates.eigenStates.eigenVectors.resize(resultsFromDavidson.eigenVectors.rows(),
                                                  resultsFromDavidson.eigenVectors.cols());
    excitedStates.eigenStates.eigenVectors.colwise() = input_.energyDifferences();
    // Loop needed in stead of rowwise needed because of a known bug in the MKL
    // assignment in Eigen/3.2.2, Bug 1527 in changelog.
    // Original version:
    // excitedStates.eigenStates.eigenVectors.array().rowwise() /=
    // resultsFromDavidson.eigenValues.array().transpose().sqrt();
    for (int row = 0; row < int(excitedStates.eigenStates.eigenVectors.rows()); ++row) {
      excitedStates.eigenStates.eigenVectors.row(row).array() /= resultsFromDavidson.eigenValues.array().transpose().sqrt();
    }
    excitedStates.eigenStates.eigenVectors = excitedStates.eigenStates.eigenVectors.array().sqrt();
    excitedStates.eigenStates.eigenVectors.array() *= resultsFromDavidson.eigenVectors.array();

    // The excitation energy is the square root of the eigenvalue
    // Temporary needed because of a known bug in the MKL assignment in Eigen/3.2.2
    // Bug 1527 in changelog.
    // Original version:
    // excitedStates.eigenStates.eigenValues = resultsFromDavidson.eigenValues.cwiseSqrt();
    Eigen::MatrixXd temp = resultsFromDavidson.eigenValues.cwiseSqrt();
    excitedStates.eigenStates.eigenValues = std::move(temp);
  }

  generateTransitionDipoleMoments(excitedStates, spinBlock);

  // generate phase of wavefunction

  Eigen::MatrixXd sign =
      Eigen::MatrixXd::Ones(excitedStates.eigenStates.eigenVectors.rows(), excitedStates.eigenStates.eigenVectors.cols());
  sign = (sign.array() > 0).select(sign, -sign);

  // Normalize excitedStates coefficient to the square coefficient sum of 1.
  excitedStates.eigenStates.eigenVectors =
      excitedStates.eigenStates.eigenVectors.array().square().matrix().cwiseProduct(sign).colwise().normalized();

  return excitedStates;
}

template<Utils::Reference restrictedness>
void TDDFTBEigenvalueSolver<restrictedness>::generateTransitionDipoleMoments(Utils::ElectronicTransitionResult& excitedStatesResults,
                                                                             Utils::SpinTransition spinBlock) const {
  // Triplets have a transition dipole moment of 0, the function transitionDipoleMoment just multiplies by 2 in the
  // restricted case so care must be taken.
  if (!dipoleMatrix_ || spinBlock == Utils::SpinTransition::Triplet) {
    excitedStatesResults.transitionDipoles = Eigen::Matrix3Xd::Zero(3, excitedStatesResults.eigenStates.eigenValues.size());
  }
  else {
    excitedStatesResults.transitionDipoles = Utils::TransitionDipoleCalculator::calculate<restrictedness>(
        *dipoleMatrix_, excitedStatesResults.eigenStates.eigenVectors, input_.excitations());
  }
}

template<>
inline auto TDDFTBEigenvalueSolver<Utils::Reference::Restricted>::generateGuess(int nConfigurations, Utils::SpinTransition spinBlock)
    -> boost::optional<Eigen::MatrixXd> {
  if (guess_) {
    if (spinBlock == Utils::SpinTransition::Singlet) {
      if (guess_->singlet.rows() != nConfigurations) {
        throw std::runtime_error("Number of configurations is not the same as the row dimention of the initial guess.");
      }
      return guess_->singlet;
    }
    else {
      if (guess_->triplet.rows() != nConfigurations) {
        throw std::runtime_error("Number of configurations is not the same as the row dimention of the initial guess.");
      }
      return guess_->triplet;
    }
  }
  return {};
}
template<>
inline auto TDDFTBEigenvalueSolver<Utils::Reference::Unrestricted>::generateGuess(int nConfigurations,
                                                                                  Utils::SpinTransition /*spinBlock*/)
    -> boost::optional<Eigen::MatrixXd> {
  if (guess_) {
    if (guess_->unrestricted.rows() != nConfigurations) {
      throw std::runtime_error("Number of configurations is not the same as the row dimention of the initial guess.");
    }
    return guess_->unrestricted;
  }
  return {};
}
} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_TDDFTBEIGENVALUESOLVER_H
