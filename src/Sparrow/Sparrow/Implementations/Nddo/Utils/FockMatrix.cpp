/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "FockMatrix.h"
#include <Sparrow/Implementations/Nddo/Utils/NDDOElectronicEnergyCalculator.h>
#include <Utils/MethodEssentials/Methods/OverlapCalculator.h>
#include <Utils/MethodEssentials/util/DensityMatrix.h>

namespace Scine {
namespace Sparrow {

namespace nddo {

FockMatrix::FockMatrix(const Utils::ElementTypeCollection& elements, const Utils::PositionCollection& positions,
                       const Utils::DensityMatrix& densityMatrix, const OneCenterIntegralContainer& oneCIntegrals,
                       const ElementParameters& elementPar, const Utils::AtomsOrbitalsIndexes& aoIndexes,
                       const Utils::OverlapCalculator& overlapCalculator, const bool& unrestrictedCalculationRunning)
  : twoCenterIntegrals_(elements, positions, elementPar),
    F1_(elements, positions, densityMatrix.restrictedMatrix(), twoCenterIntegrals_, elementPar, aoIndexes),
    F2_(elements, densityMatrix, oneCIntegrals, twoCenterIntegrals_, elementPar, aoIndexes),
    overlapCalculator_(overlapCalculator),
    unrestrictedCalculationRunning_(unrestrictedCalculationRunning) {
  electronicEnergyCalculator_ =
      std::make_unique<NDDOElectronicEnergyCalculator>(densityMatrix, *this, unrestrictedCalculationRunning_);
}

void FockMatrix::initialize() {
  twoCenterIntegrals_.initialize();
  F1_.initialize();
  F2_.initialize();
}

void FockMatrix::calculateDensityIndependentPart(Utils::derivOrder order) {
  twoCenterIntegrals_.update(order);
  F1_.calculate(overlapCalculator_.getOverlap()); // NEEDS TO BE AFTER twoCenterIntegrals update!
}

void FockMatrix::calculateDensityDependentPart(Utils::derivOrder /*order*/) {
  F2_.calculate(unrestrictedCalculationRunning_);
}

void FockMatrix::addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::first>& derivatives) const {
  addDerivativesImpl<Utils::derivativeType::first>(derivatives);
}

void FockMatrix::addDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives) const {
  addDerivativesImpl<Utils::derivativeType::second_atomic>(derivatives);
}

void FockMatrix::addDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_full>& derivatives) const {
  addDerivativesImpl<Utils::derivativeType::second_full>(derivatives);
}

template<Utils::derivativeType O>
void FockMatrix::addDerivativesImpl(Utils::AutomaticDifferentiation::DerivativeContainerType<O>& derivatives) const {
  F1_.addDerivatives<O>(derivatives, overlapCalculator_.getOverlap());
  F2_.addDerivatives<O>(derivatives);
}

const OneElectronMatrix& FockMatrix::getOneElectronMatrix() const {
  return F1_;
}

const TwoElectronMatrix& FockMatrix::getTwoElectronMatrix() const {
  return F2_;
}

Utils::SpinAdaptedMatrix FockMatrix::getMatrix() const {
  Utils::SpinAdaptedMatrix fock;
  if (!unrestrictedCalculationRunning_)
    fock.setRestrictedMatrix(F1_.getMatrix() + F2_.getMatrix());
  else {
    fock.setAlphaMatrix(F1_.getMatrix() + F2_.getAlpha());
    fock.setBetaMatrix(F1_.getMatrix() + F2_.getBeta());
  }
  return fock;
}

double FockMatrix::calculateElectronicEnergy() const {
  return electronicEnergyCalculator_->calculateElectronicEnergy();
}

void FockMatrix::finalize(Utils::derivOrder order) {
  // Recalculate the Fock matrix: make it consistent with the obtained density matrix
  calculateDensityDependentPart(order);
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
