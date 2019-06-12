/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ScfFock.h"
#include "ZeroOrderMatricesCalculator.h"
#include <Utils/MethodEssentials/util/DensityMatrix.h>
#include <Utils/MethodEssentials/util/LcaoUtil/LcaoUtil.h>

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

void ScfFock::calculateDensityDependentPart(Utils::derivOrder /*order*/) {
  populationAnalysis();
  if (unrestrictedCalculationRunning_) {
    spinDFTB.spinPopulationAnalysis(densityMatrix_.alphaMatrix(), densityMatrix_.betaMatrix(), overlapMatrix_);
    spinDFTB.calculateSpinContribution();
  }
  completeH();
}

void ScfFock::calculateDensityIndependentPart(Utils::derivOrder order) {
  zeroOrderMatricesCalculator_.calculateFockMatrix(order);
  H0_ = zeroOrderMatricesCalculator_.getZeroOrderHamiltonian().getMatrixXd();
  constructG(order);
}

Utils::SpinAdaptedMatrix ScfFock::getMatrix() const {
  Eigen::MatrixXd sum = zeroOrderMatricesCalculator_.getZeroOrderHamiltonian().getMatrixXd() + correctionToFock;
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
  Utils::LcaoUtil::calculateMullikenAtomicCharges(atomicCharges_, coreCharges_, densityMatrix_, overlapMatrix_, aoIndexes_);
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
