/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "TDDFTBSigmaVectorEvaluator.h"
#include <Sparrow/Implementations/Dftb/TimeDependent/LinearResponse/BasisPruner.h>
#include <Sparrow/Implementations/Exceptions.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <iostream>

namespace Scine {
namespace Sparrow {

namespace detail {

inline void negativeBetaTransitions(Eigen::Ref<Eigen::VectorXd> toTransform,
                                    const Eigen::Ref<const Eigen::Matrix<bool, -1, 1>>& isBeta) {
  assert(toTransform.size() == isBeta.size());
  toTransform = (isBeta).select(-toTransform, toTransform);
}
} // namespace detail

template<Utils::Reference restrictedness>
TDDFTBSigmaVectorEvaluator<restrictedness>::TDDFTBSigmaVectorEvaluator(std::shared_ptr<Eigen::MatrixXd> gammaMatrix,
                                                                       std::shared_ptr<Eigen::VectorXd> spinConstants,
                                                                       const OrderedInput<restrictedness>& orderedInput,
                                                                       Utils::SpinTransition spinBlock, TDDFTBType type)
  : input_(orderedInput),
    gammaMatrix_(std::move(gammaMatrix)),
    spinConstantsVector_(std::move(spinConstants)),
    isTDA_(type == TDDFTBType::TDA),
    spinBlock_(spinBlock) {
  check();
  calculateAtomicEnergyWeightedTransitionCharges(input_.transitionCharges());
}

template<Utils::Reference restrictedness>
const Eigen::MatrixXd& TDDFTBSigmaVectorEvaluator<restrictedness>::evaluate(const Eigen::MatrixXd& guessVectors) const {
  assert(guessVectors.rows() == input_.energyDifferences().size());

  const int dimCol = guessVectors.cols();
  const int dimRow = guessVectors.rows();

  const int alreadyComputedSigmaVectors = currentSigmaMatrix_.cols();
  const int vectorsToCompute = dimCol - alreadyComputedSigmaVectors;

  Eigen::MatrixXd sigmaMatrix(dimRow, vectorsToCompute);
  // X_{BI} = \sum_{jb} h_{jb, B} * T_{jb, I}
  // T: guess vector
  // h: atomic energy weighted transition charges
  // Y_{AI} = \sum_{B} \gamma_{AB} * X_{BI}
  // R_{ia,I} = \Delta_{ia}^2 * T_{ia, I} + 4 * \sum_{A} h_{ia, A} * Y_{AI}
  // R: sigma vector

#pragma omp parallel for schedule(dynamic) default(none) shared(sigmaMatrix, guessVectors) \
    firstprivate(vectorsToCompute, alreadyComputedSigmaVectors)
  for (int col = 0; col < vectorsToCompute; ++col) {
    int const colIndex = col + alreadyComputedSigmaVectors;
    const Eigen::VectorXd& guessVector = guessVectors.col(colIndex);

    // Necessary to have a temporary here because with some
    // setup of linear algebra libraries not having this caused race conditions.
    Eigen::VectorXd sigmaMatrixColumn =
        isTDA_ ? Eigen::VectorXd(input_.energyDifferences().cwiseProduct(guessVector))
               : Eigen::VectorXd(input_.energyDifferences().cwiseAbs2().cwiseProduct(guessVector));

    sigmaMatrixColumn += calculateAtomicContraction(calculateYAI(calculateXBI(guessVector)));

    fillAdditionalSigmaMatrixTerms(sigmaMatrixColumn, guessVector);
    sigmaMatrix.col(col) = sigmaMatrixColumn;
  }

  currentSigmaMatrix_.conservativeResize(dimRow, dimCol);
  currentSigmaMatrix_.rightCols(vectorsToCompute) = sigmaMatrix;
  return currentSigmaMatrix_;
}

template<Utils::Reference restrictedness>
void TDDFTBSigmaVectorEvaluator<restrictedness>::calculateAtomicEnergyWeightedTransitionCharges(const Eigen::MatrixXd& transitionCharges) {
  assert(transitionCharges.rows() == input_.energyDifferences().rows());
  if (isTDA_) {
    energyWeightedAtomicTransitionCharges_ = transitionCharges;
  }
  else {
    energyWeightedAtomicTransitionCharges_ =
        transitionCharges.array().colwise() * input_.energyDifferences().cwiseSqrt().array();
  }
}

template<Utils::Reference restrictedness>
Eigen::MatrixXd TDDFTBSigmaVectorEvaluator<restrictedness>::calculateXBI(const Eigen::VectorXd& guessVector) const {
  return energyWeightedAtomicTransitionCharges_.transpose() * guessVector;
}

template<>
Eigen::MatrixXd TDDFTBSigmaVectorEvaluator<Utils::Reference::Restricted>::calculateYAI(const Eigen::VectorXd& XBI) const {
  assert(spinBlock_ == Utils::SpinTransition::Singlet || spinBlock_ == Utils::SpinTransition::Triplet);
  if (spinBlock_ == Utils::SpinTransition::Singlet) {
    return gammaMatrix_->selfadjointView<Eigen::Lower>() * XBI;
  }
  else { // Utils::SpinTransition::Triplet
    return spinConstantsVector_->cwiseProduct(XBI);
  }
}

template<>
Eigen::MatrixXd TDDFTBSigmaVectorEvaluator<Utils::Reference::Unrestricted>::calculateYAI(const Eigen::VectorXd& XBI) const {
  return (gammaMatrix_->selfadjointView<Eigen::Lower>()) * XBI;
}

template<Utils::Reference restrictedness>
Eigen::VectorXd TDDFTBSigmaVectorEvaluator<restrictedness>::calculateAtomicContraction(const Eigen::MatrixXd& YAIMatrix) const {
  return factor() * energyWeightedAtomicTransitionCharges_ * YAIMatrix;
}

template<Utils::Reference restrictedness>
void TDDFTBSigmaVectorEvaluator<restrictedness>::collapsed(int /*newSubspaceDimension*/) {
  currentSigmaMatrix_ = Eigen::MatrixXd(0, 0);
}

template<>
template<typename Derived, typename OtherDerived>
inline void TDDFTBSigmaVectorEvaluator<Utils::Reference::Restricted>::fillAdditionalSigmaMatrixTerms(
    const Eigen::MatrixBase<Derived>& /*sigmaBlock*/, const Eigen::MatrixBase<OtherDerived>& /*guessVector*/) const {
}

template<>
template<typename Derived, typename OtherDerived>
inline void TDDFTBSigmaVectorEvaluator<Utils::Reference::Unrestricted>::fillAdditionalSigmaMatrixTerms(
    const Eigen::MatrixBase<Derived>& sigmaBlock, const Eigen::MatrixBase<OtherDerived>& guessVector) const {
  assert(guessVector.size() == input_.isBeta().size());
  Eigen::VectorXd modifiedGuessVector = guessVector;

  detail::negativeBetaTransitions(modifiedGuessVector, input_.isBeta());

  for (int atom = 0; atom < spinConstantsVector_->size(); ++atom) {
    double spinDependentTransitionCont = energyWeightedAtomicTransitionCharges_.col(atom).transpose() * modifiedGuessVector;
    spinDependentTransitionCont *= 2.0 * (*spinConstantsVector_)(atom);

    Eigen::VectorXd spinDependentTransition = spinDependentTransitionCont * energyWeightedAtomicTransitionCharges_.col(atom);

    detail::negativeBetaTransitions(spinDependentTransition, input_.isBeta());

    const_cast<Eigen::MatrixBase<Derived>&>(sigmaBlock) += spinDependentTransition;
  }
}

template<>
void TDDFTBSigmaVectorEvaluator<Utils::Reference::Restricted>::check() const {
  assert(spinBlock_ == Utils::SpinTransition::Singlet || spinBlock_ == Utils::SpinTransition::Triplet);
  if (spinBlock_ == Utils::SpinTransition::Triplet && !spinConstantsVector_) {
    throw SpinConstantsNotAvailableException();
  }
}

template<>
void TDDFTBSigmaVectorEvaluator<Utils::Reference::Unrestricted>::check() const {
  if (!spinConstantsVector_) {
    throw SpinConstantsNotAvailableException();
  }
}

template class TDDFTBSigmaVectorEvaluator<Utils::Reference::Unrestricted>;
template class TDDFTBSigmaVectorEvaluator<Utils::Reference::Restricted>;
} // namespace Sparrow
} // namespace Scine
