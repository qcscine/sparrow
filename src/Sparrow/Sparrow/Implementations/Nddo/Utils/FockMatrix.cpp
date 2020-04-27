/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "FockMatrix.h"
#include <Sparrow/Implementations/Nddo/Utils/NDDOElectronicEnergyCalculator.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/Scf/MethodInterfaces/AdditiveElectronicContribution.h>
#include <Utils/Scf/MethodInterfaces/OverlapCalculator.h>

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
  for (auto& contribution : densityIndependentContributions_) {
    if (contribution->isValid())
      contribution->calculate({}, order);
  }
}

void FockMatrix::calculateDensityDependentPart(Utils::derivOrder order) {
  F2_.calculate(unrestrictedCalculationRunning_);

  for (auto& contribution : densityDependentContributions_) {
    if (contribution->isValid())
      contribution->calculate({}, order);
  }
}

void FockMatrix::addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::first>& derivatives) const {
  addDerivativesImpl<Utils::derivativeType::first>(derivatives);
  for (auto& contribution : densityIndependentContributions_) {
    if (contribution->isValid())
      contribution->addDerivatives(derivatives);
  }
  for (auto& contribution : densityDependentContributions_) {
    if (contribution->isValid())
      contribution->addDerivatives(derivatives);
  }
}

void FockMatrix::addDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives) const {
  addDerivativesImpl<Utils::derivativeType::second_atomic>(derivatives);
  for (auto& contribution : densityIndependentContributions_) {
    if (contribution->isValid())
      contribution->addDerivatives(derivatives);
  }
  for (auto& contribution : densityDependentContributions_) {
    if (contribution->isValid())
      contribution->addDerivatives(derivatives);
  }
}

void FockMatrix::addDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_full>& derivatives) const {
  addDerivativesImpl<Utils::derivativeType::second_full>(derivatives);
  for (auto& contribution : densityIndependentContributions_) {
    if (contribution->isValid())
      contribution->addDerivatives(derivatives);
  }
  for (auto& contribution : densityDependentContributions_) {
    if (contribution->isValid())
      contribution->addDerivatives(derivatives);
  }
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
  if (!unrestrictedCalculationRunning_) {
    Eigen::MatrixXd restrictedFock = F1_.getMatrix() + F2_.getMatrix();
    for (const auto& contribution : densityDependentContributions_) {
      if (contribution->isValid() && contribution->hasMatrixContribution())
        restrictedFock += contribution->getElectronicContribution().restrictedMatrix();
    }
    for (const auto& contribution : densityIndependentContributions_) {
      if (contribution->isValid() && contribution->hasMatrixContribution()) {
        restrictedFock += contribution->getElectronicContribution().restrictedMatrix();
      }
    }
    fock.setRestrictedMatrix(std::move(restrictedFock));
  }
  else {
    Eigen::MatrixXd unrestrictedFock = F1_.getMatrix();
    for (const auto& contribution : densityDependentContributions_) {
      if (contribution->isValid() && contribution->hasMatrixContribution()) {
        unrestrictedFock += contribution->getElectronicContribution().alphaMatrix() +
                            contribution->getElectronicContribution().betaMatrix();
      }
    }
    for (const auto& contribution : densityIndependentContributions_) {
      if (contribution->isValid() && contribution->hasMatrixContribution()) {
        unrestrictedFock += contribution->getElectronicContribution().alphaMatrix() +
                            contribution->getElectronicContribution().betaMatrix();
      }
    }
    fock.setAlphaMatrix(unrestrictedFock + F2_.getAlpha());
    fock.setBetaMatrix(unrestrictedFock + F2_.getBeta());
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

void FockMatrix::addDensityIndependentElectronicContribution(std::shared_ptr<Utils::AdditiveElectronicContribution> contribution) {
  densityIndependentContributions_.emplace_back(std::move(contribution));
}

void FockMatrix::addDensityDependentElectronicContribution(std::shared_ptr<Utils::AdditiveElectronicContribution> contribution) {
  densityDependentContributions_.emplace_back(std::move(contribution));
}

const std::vector<std::shared_ptr<Utils::AdditiveElectronicContribution>>& FockMatrix::getDensityDependentContributions() const {
  return densityDependentContributions_;
}

const std::vector<std::shared_ptr<Utils::AdditiveElectronicContribution>>& FockMatrix::getDensityIndependentContributions() const {
  return densityIndependentContributions_;
}

void FockMatrix::clearElectronicContributions() {
  densityIndependentContributions_.clear();
  densityIndependentContributions_.clear();
}

void FockMatrix::eraseElectronicContribution(std::shared_ptr<Utils::AdditiveElectronicContribution> contribution) {
  auto it = std::find(densityIndependentContributions_.begin(), densityIndependentContributions_.end(), contribution);
  if (it != densityIndependentContributions_.end()) {
    densityIndependentContributions_.erase(it);
  }
  it = std::find(densityDependentContributions_.begin(), densityDependentContributions_.end(), contribution);
  if (it != densityDependentContributions_.end()) {
    densityDependentContributions_.erase(it);
  }
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
