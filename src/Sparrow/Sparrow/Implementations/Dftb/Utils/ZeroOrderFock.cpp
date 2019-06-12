/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ZeroOrderFock.h"
#include "ZeroOrderMatricesCalculator.h"
#include <Utils/MethodEssentials/util/SingleParticleEnergies.h>
#include <Utils/MethodEssentials/util/SpinAdaptedMatrix.h>

namespace Scine {
namespace Sparrow {

namespace dftb {

ZeroOrderFock::ZeroOrderFock(ZeroOrderMatricesCalculator& matricesCalculator,
                             const Utils::SingleParticleEnergies& singleParticleEnergies,
                             const Eigen::MatrixXd& energyWeightedDensityMatrix, const int& nElectrons)
  : matricesCalculator_(matricesCalculator),
    singleParticleEnergies_(singleParticleEnergies),
    energyWeightedDensityMatrix_(energyWeightedDensityMatrix),
    nElectrons_(nElectrons) {
}

void ZeroOrderFock::initialize() {
  matricesCalculator_.initializeFockCalculator();
}

void ZeroOrderFock::calculateDensityIndependentPart(Utils::derivOrder order) {
  matricesCalculator_.calculateFockMatrix(order);
}

void ZeroOrderFock::calculateDensityDependentPart(Utils::derivOrder /*order*/) {
  // Nothing: only density-independent.
}

Utils::SpinAdaptedMatrix ZeroOrderFock::getMatrix() const {
  Utils::SpinAdaptedMatrix fock;
  fock.setRestrictedMatrix(matricesCalculator_.getZeroOrderHamiltonian().getMatrixXd());
  return fock;
}

void ZeroOrderFock::addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::first>& derivatives) const {
  matricesCalculator_.addDerivatives(derivatives, energyWeightedDensityMatrix_);
}

void ZeroOrderFock::addDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives) const {
  matricesCalculator_.addDerivatives(derivatives, energyWeightedDensityMatrix_);
}

void ZeroOrderFock::addDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_full>& derivatives) const {
  matricesCalculator_.addDerivatives(derivatives, energyWeightedDensityMatrix_);
}

double ZeroOrderFock::calculateElectronicEnergy() const {
  double elEnergy = 0.0;
  for (int i = 0; i < nElectrons_ / 2; ++i)
    elEnergy += singleParticleEnergies_.getRestrictedEnergies()[i] * 2;

  // If there is a singly occupied orbital, add its energy
  if (nElectrons_ % 2 != 0)
    elEnergy += singleParticleEnergies_.getRestrictedEnergies()[nElectrons_ / 2];

  return elEnergy;
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
