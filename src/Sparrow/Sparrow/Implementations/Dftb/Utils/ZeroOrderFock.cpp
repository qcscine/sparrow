/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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

void ZeroOrderFock::calculateDensityIndependentPart(Utils::DerivativeOrder order) {
  matricesCalculator_.calculateFockMatrix(order);
  for (auto& contribution : densityIndependentContributions_) {
    contribution->calculate({}, order);
  }
}

void ZeroOrderFock::calculateDensityDependentPart(Utils::DerivativeOrder /*order*/) {
  // Nothing: only density-independent.
}

Utils::SpinAdaptedMatrix ZeroOrderFock::getMatrix() const {
  Utils::SpinAdaptedMatrix fock;
  fock.setRestrictedMatrix(matricesCalculator_.getZeroOrderHamiltonian().getMatrixXd());
  for (const auto& contribution : densityIndependentContributions_) {
    if (contribution->isValid() && contribution->hasMatrixContribution()) {
      fock.restrictedMatrix() += contribution->getElectronicContribution().restrictedMatrix();
    }
  }
  return fock;
}

void ZeroOrderFock::addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::Derivative::First>& derivatives) const {
  matricesCalculator_.addDerivatives(derivatives, energyWeightedDensityMatrix_);
  for (const auto& contribution : densityIndependentContributions_) {
    if (contribution->isValid()) {
      contribution->addDerivatives(derivatives);
    }
  }
}

void ZeroOrderFock::addDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::Derivative::SecondAtomic>& derivatives) const {
  matricesCalculator_.addDerivatives(derivatives, energyWeightedDensityMatrix_);
  for (const auto& contribution : densityIndependentContributions_) {
    if (contribution->isValid()) {
      contribution->addDerivatives(derivatives);
    }
  }
}

void ZeroOrderFock::addDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::Derivative::SecondFull>& derivatives) const {
  matricesCalculator_.addDerivatives(derivatives, energyWeightedDensityMatrix_);
  for (const auto& contribution : densityIndependentContributions_) {
    if (contribution->isValid()) {
      contribution->addDerivatives(derivatives);
    }
  }
}

double ZeroOrderFock::calculateElectronicEnergy() const {
  const auto& restrictedEnergies = singleParticleEnergies_.getRestrictedEnergies();

  double elEnergy = 0.0;
  for (int i = 0; i < nElectrons_ / 2; ++i) {
    elEnergy += restrictedEnergies[i] * 2;
  }

  // If there is a singly occupied orbital, add its energy
  if (nElectrons_ % 2 != 0) {
    elEnergy += restrictedEnergies[nElectrons_ / 2];
  }

  for (const auto& contribution : densityIndependentContributions_) {
    elEnergy += contribution->getElectronicEnergyContribution();
  }

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
