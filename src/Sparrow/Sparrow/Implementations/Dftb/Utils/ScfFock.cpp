/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ScfFock.h"
#include "ZeroOrderMatricesCalculator.h"
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/Scf/LcaoUtils/LcaoUtils.h>
#include <Utils/Scf/MethodInterfaces/AdditiveElectronicContribution.h>

namespace Scine {
namespace Sparrow {

namespace dftb {

ScfFock::ScfFock(ZeroOrderMatricesCalculator& matricesCalculator, const Utils::ElementTypeCollection& elements,
                 const Utils::PositionCollection& positions, const DFTBCommon::AtomicParameterContainer& atomicPar,
                 const DFTBCommon::DiatomicParameterContainer& diatomicPar, const Utils::DensityMatrix& densityMatrix,
                 const Eigen::MatrixXd& energyWeightedDensityMatrix, std::vector<double>& atomicCharges,
                 const std::vector<double>& coreCharges, const Utils::AtomsOrbitalsIndexes& aoIndexes,
                 const Eigen::MatrixXd& overlapMatrix, const bool& unrestrictedCalculationRunning)
  : zeroOrderMatricesCalculator_(matricesCalculator),
    elements_(elements),
    positions_(positions),
    atomicPar_(atomicPar),
    diatomicPar_(diatomicPar),
    densityMatrix_(densityMatrix),
    energyWeightedDensityMatrix_(energyWeightedDensityMatrix),
    atomicCharges_(atomicCharges),
    coreCharges_(coreCharges),
    aoIndexes_(aoIndexes),
    overlapMatrix_(overlapMatrix),
    unrestrictedCalculationRunning_(unrestrictedCalculationRunning),
    spinDFTB(elements, atomicPar) {
}

void ScfFock::initialize() {
  auto numberOrbitals = aoIndexes_.getNAtomicOrbitals();

  zeroOrderMatricesCalculator_.initializeFockCalculator();
  HXoverS_ = Eigen::MatrixXd::Zero(numberOrbitals, numberOrbitals);
  correctionToFock = Eigen::MatrixXd::Zero(numberOrbitals, numberOrbitals);

  spinDFTB.initialize(getNumberAtoms(), numberOrbitals, aoIndexes_);
}

void ScfFock::calculateDensityDependentPart(Utils::derivOrder order) {
  populationAnalysis();
  if (unrestrictedCalculationRunning_) {
    spinDFTB.spinPopulationAnalysis(densityMatrix_.alphaMatrix(), densityMatrix_.betaMatrix(), overlapMatrix_);
    spinDFTB.calculateSpinContribution();
  }
  completeH();
  for (auto& contribution : densityDependentContributions_) {
    contribution->calculate(densityMatrix_, order);
  }
  for (auto& contribution : densityIndependentContributions_) {
    contribution->calculate(densityMatrix_, order);
  }
}

void ScfFock::calculateDensityIndependentPart(Utils::derivOrder order) {
  zeroOrderMatricesCalculator_.calculateFockMatrix(order);
  H0_ = zeroOrderMatricesCalculator_.getZeroOrderHamiltonian().getMatrixXd();
  constructG(order);
  for (auto& contribution : densityDependentContributions_) {
    contribution->calculate(densityMatrix_, order);
  }
  for (auto& contribution : densityIndependentContributions_) {
    contribution->calculate(densityMatrix_, order);
  }
}

Utils::SpinAdaptedMatrix ScfFock::getMatrix() const {
  Eigen::MatrixXd sum = zeroOrderMatricesCalculator_.getZeroOrderHamiltonian().getMatrixXd() + correctionToFock;
  for (auto& contribution : densityDependentContributions_) {
    sum += contribution->getElectronicContribution().restrictedMatrix();
  }
  for (auto& contribution : densityIndependentContributions_) {
    contribution->getElectronicContribution().restrictedMatrix();
  }

  Utils::SpinAdaptedMatrix fock;
  fock.setRestrictedMatrix(std::move(sum));
  if (unrestrictedCalculationRunning_)
    spinDFTB.constructSpinHamiltonians(fock, overlapMatrix_);
  return fock;
}

void ScfFock::finalize(Utils::derivOrder /*order*/) {
  // Repeat the population analysis to make sure that the correct spin energy is employed
  // populationAnalysis(); (not needed anymore; is already done in ScfMethod::finalizeCalculation.
  if (unrestrictedCalculationRunning_) {
    spinDFTB.spinPopulationAnalysis(densityMatrix_.alphaMatrix(), densityMatrix_.betaMatrix(), overlapMatrix_);
  }
}

void ScfFock::populationAnalysis() {
  Utils::LcaoUtils::calculateMullikenAtomicCharges(atomicCharges_, coreCharges_, densityMatrix_, overlapMatrix_, aoIndexes_);
}

void ScfFock::addDensityDependentElectronicContribution(std::shared_ptr<Utils::AdditiveElectronicContribution> contribution) {
  densityDependentContributions_.emplace_back(std::move(contribution));
}

void ScfFock::addDensityIndependentElectronicContribution(std::shared_ptr<Utils::AdditiveElectronicContribution> contribution) {
  densityIndependentContributions_.emplace_back(std::move(contribution));
}

void ScfFock::addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::first>& derivatives) const {
  for (auto& contribution : densityIndependentContributions_) {
    if (contribution->isValid())
      contribution->addDerivatives(derivatives);
  }
  for (auto& contribution : densityDependentContributions_) {
    if (contribution->isValid())
      contribution->addDerivatives(derivatives);
  }
}

void ScfFock::addDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives) const {
  for (auto& contribution : densityIndependentContributions_) {
    if (contribution->isValid())
      contribution->addDerivatives(derivatives);
  }
  for (auto& contribution : densityDependentContributions_) {
    if (contribution->isValid())
      contribution->addDerivatives(derivatives);
  }
}

void ScfFock::addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_full>& derivatives) const {
  for (auto& contribution : densityIndependentContributions_) {
    if (contribution->isValid())
      contribution->addDerivatives(derivatives);
  }
  for (auto& contribution : densityDependentContributions_) {
    if (contribution->isValid())
      contribution->addDerivatives(derivatives);
  }
}
} // namespace dftb
} // namespace Sparrow
} // namespace Scine
