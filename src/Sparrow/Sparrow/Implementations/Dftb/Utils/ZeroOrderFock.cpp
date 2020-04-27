/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ZeroOrderFock.h"
#include "ZeroOrderMatricesCalculator.h"
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/SingleParticleEnergies.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <Utils/Scf/MethodInterfaces/AdditiveElectronicContribution.h>

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
  for (auto& contribution : densityIndependentContributions_) {
    contribution->calculate({}, order);
  }
}

void ZeroOrderFock::calculateDensityDependentPart(Utils::derivOrder /*order*/) {
  // Nothing: only density-independent.
}

Utils::SpinAdaptedMatrix ZeroOrderFock::getMatrix() const {
  Utils::SpinAdaptedMatrix fock;
  fock.setRestrictedMatrix(matricesCalculator_.getZeroOrderHamiltonian().getMatrixXd());
  for (auto& contribution : densityIndependentContributions_) {
    if (contribution->isValid() && contribution->hasMatrixContribution())
      fock.restrictedMatrix() += contribution->getElectronicContribution().restrictedMatrix();
  }
  return fock;
}

void ZeroOrderFock::addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::first>& derivatives) const {
  matricesCalculator_.addDerivatives(derivatives, energyWeightedDensityMatrix_);
  for (auto& contribution : densityIndependentContributions_) {
    if (contribution->isValid())
      contribution->addDerivatives(derivatives);
  }
}

void ZeroOrderFock::addDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives) const {
  matricesCalculator_.addDerivatives(derivatives, energyWeightedDensityMatrix_);
  for (auto& contribution : densityIndependentContributions_) {
    if (contribution->isValid())
      contribution->addDerivatives(derivatives);
  }
}

void ZeroOrderFock::addDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_full>& derivatives) const {
  matricesCalculator_.addDerivatives(derivatives, energyWeightedDensityMatrix_);
  for (auto& contribution : densityIndependentContributions_) {
    if (contribution->isValid())
      contribution->addDerivatives(derivatives);
  }
}

double ZeroOrderFock::calculateElectronicEnergy() const {
  double elEnergy = 0.0;
  for (int i = 0; i < nElectrons_ / 2; ++i)
    elEnergy += singleParticleEnergies_.getRestrictedEnergies()[i] * 2;

  // If there is a singly occupied orbital, add its energy
  if (nElectrons_ % 2 != 0)
    elEnergy += singleParticleEnergies_.getRestrictedEnergies()[nElectrons_ / 2];

  for (auto& contribution : densityIndependentContributions_)
    elEnergy += contribution->getElectronicEnergyContribution();

  return elEnergy;
}

void ZeroOrderFock::addDensityDependentElectronicContribution(std::shared_ptr<Utils::AdditiveElectronicContribution> contribution) {
  densityIndependentContributions_.emplace_back(std::move(contribution));
}

void ZeroOrderFock::addDensityIndependentElectronicContribution(std::shared_ptr<Utils::AdditiveElectronicContribution> contribution) {
  densityIndependentContributions_.emplace_back(std::move(contribution));
}
} // namespace dftb
} // namespace Sparrow
} // namespace Scine
