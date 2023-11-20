/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "DiagonalPreconditionerEvaluator.h"
#include "TimeDependentUtils.h"

namespace Scine {
namespace Sparrow {

DiagonalPreconditionerEvaluator::DiagonalPreconditionerEvaluator(const Eigen::VectorXd& energyDifferenceVector) {
  energyDifferences_ = energyDifferenceVector;
}

DiagonalPreconditionerEvaluator::DiagonalPreconditionerEvaluator(
    const Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::VectorXd>& energyDifferenceVector) {
  int nAlpha = energyDifferenceVector.alpha.size();
  int nBeta = energyDifferenceVector.beta.size();
  energyDifferences_.resize(nAlpha + nBeta);
  energyDifferences_.head(nAlpha) = energyDifferenceVector.alpha;
  energyDifferences_.tail(nBeta) = energyDifferenceVector.beta;
}

DiagonalPreconditionerEvaluator::DiagonalPreconditionerEvaluator(
    const Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd>& energyDifferenceVector) {
  energyDifferences_ = energyDifferenceVector.restricted;
}

DiagonalPreconditionerEvaluator::DiagonalPreconditionerEvaluator(
    const Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::VectorXd>& energyDifferenceVector, OrderTag) {
  int nAlpha = energyDifferenceVector.alpha.size();
  int nBeta = energyDifferenceVector.beta.size();
  Eigen::VectorXd vectorInStandardOrder(nAlpha + nBeta);
  vectorInStandardOrder.head(nAlpha) = energyDifferenceVector.alpha;
  vectorInStandardOrder.tail(nBeta) = energyDifferenceVector.beta;
  TimeDependentUtils::transformOrder(vectorInStandardOrder, energyDifferences_,
                                     TimeDependentUtils::generateEnergyOrderMap(energyDifferenceVector),
                                     TimeDependentUtils::Direction::To);
}

DiagonalPreconditionerEvaluator::DiagonalPreconditionerEvaluator(
    const Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd>& energyDifferenceVector, OrderTag) {
  TimeDependentUtils::transformOrder(energyDifferenceVector.restricted, energyDifferences_,
                                     TimeDependentUtils::generateEnergyOrderMap(energyDifferenceVector),
                                     TimeDependentUtils::Direction::To);
}

Eigen::VectorXd DiagonalPreconditionerEvaluator::evaluate(const Eigen::VectorXd& vectorToPrecondition, double eigenvalue) const {
  assert(energyDifferences_.size() == vectorToPrecondition.size());
  Eigen::VectorXd difference = eigenvalue - energyDifferences_.array();
  Eigen::VectorXd preconditionedVector = vectorToPrecondition;
  for (int i = 0; i < vectorToPrecondition.size(); ++i) {
    // From numerical experiment. Seems to avoid lack of convergence if eigenvalue is too
    // close to diagonal element.
    if (std::abs(difference(i)) < 1e-3) {
      preconditionedVector(i) = vectorToPrecondition(i);
    }
    else {
      preconditionedVector(i) /= difference(i);
    }
  }
  return preconditionedVector;
}

} // namespace Sparrow
} // namespace Scine
