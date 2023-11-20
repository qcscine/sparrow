/**
 * @file GuessPropagator.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "GuessPropagator.h"
#include <Eigen/Dense>

namespace Scine {
namespace Sparrow {
namespace RealTimeSpectroscopy {

auto GuessPropagator::calculateGuessAtNewPosition(const Utils::PositionCollection& newPositions)
    -> std::shared_ptr<LinearResponseCalculator::GuessSpecifier> {
  LinearResponseCalculator::GuessSpecifier guess;
  if (excitedStatesHistory_.empty()) {
    return std::shared_ptr<LinearResponseCalculator::GuessSpecifier>();
  }
  Eigen::MatrixXd& excitedStatesGuess = symmetry_ == SpectralType::Unrestricted ? guess.unrestricted : guess.singlet;

  // Generate linear system to solve
  if (excitedStatesHistory_.size() == 1) {
    excitedStatesGuess = excitedStatesHistory_.back();
  }
  else if (!positionsHistory_.empty()) {
    Eigen::MatrixXd lhsMatrix = calculateOverlapMatrix(newPositions);
    Eigen::VectorXd coefficientVector = lhsMatrix.colPivHouseholderQr().solve(calculateRhsVector());
    excitedStatesGuess = Eigen::MatrixXd::Zero(excitedStatesHistory_.front().rows(), excitedStatesHistory_.front().cols());
    int index = 0;
    for (const Eigen::MatrixXd& excitedState : excitedStatesHistory_)
      excitedStatesGuess += excitedState * coefficientVector(index++);
  }

  return std::make_shared<LinearResponseCalculator::GuessSpecifier>(std::move(guess));
}

void GuessPropagator::record(const LinearResponseCalculator::GuessSpecifier& excitedState) {
  if (excitedState.singlet.size() != 0) {
    if (symmetry_ == SpectralType::Unrestricted) {
      reset();
    }
    excitedStatesHistory_.push_front(excitedState.singlet);
    symmetry_ = SpectralType::Singlet;
  }
  else if (excitedState.unrestricted.size() != 0) {
    if (symmetry_ == SpectralType::Singlet) {
      reset();
    }
    excitedStatesHistory_.push_front(excitedState.unrestricted);
    symmetry_ = SpectralType::Unrestricted;
  }
  else {
    throw std::runtime_error("No singlet or unrestricted results.");
  }
  if (static_cast<int>(excitedStatesHistory_.size()) > maxHistoryDimension_) {
    excitedStatesHistory_.pop_back();
  }
}

void GuessPropagator::record(const Utils::PositionCollection& newPositions) {
  positionsHistory_.push_front(newPositions);
  if (static_cast<int>(positionsHistory_.size()) > maxHistoryDimension_) {
    positionsHistory_.pop_back();
  }
}

void GuessPropagator::reset() {
  excitedStatesHistory_.clear();
  positionsHistory_.clear();
}

Eigen::MatrixXd GuessPropagator::calculateErrorVectors(const Utils::PositionCollection& newPosition) const {
  Eigen::MatrixXd errorMatrixBlock(positionsHistory_.front().size(), positionsHistory_.size());
  int index = 0;
  for (const Utils::PositionCollection& position : positionsHistory_) {
    Utils::PositionCollection difference = position - newPosition;
    errorMatrixBlock.col(index++) = Eigen::Map<const Eigen::VectorXd>(difference.data(), difference.size());
  }
  return errorMatrixBlock;
}

Eigen::MatrixXd GuessPropagator::calculateOverlapMatrix(const Utils::PositionCollection& positions) const {
  Eigen::Index dimension = positionsHistory_.size() + 1;
  Eigen::MatrixXd errorOverlapMatrix = Eigen::MatrixXd::Constant(dimension, dimension, -1.0);
  errorOverlapMatrix(0, 0) = 0;
  // Create convert deque to eigen matrix
  Eigen::MatrixXd errorVectors = calculateErrorVectors(positions);
  // Calculate overlap
  errorOverlapMatrix.block(1, 1, dimension - 1, dimension - 1) = errorVectors.transpose() * errorVectors;

  return errorOverlapMatrix;
}

Eigen::VectorXd GuessPropagator::calculateRhsVector() const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(positionsHistory_.size() + 1);
  result(0) = -1;
  return result;
}

void GuessPropagator::setDimension(int diisDimension) {
  maxHistoryDimension_ = diisDimension;
}
} // namespace RealTimeSpectroscopy
} // namespace Sparrow
} // namespace Scine
