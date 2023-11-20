/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "NDDOElectronicEnergyCalculator.h"
#include "FockMatrix.h"
#include "OneElectronMatrix.h"
#include "TwoElectronMatrix.h"
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/Scf/MethodInterfaces/AdditiveElectronicContribution.h>

namespace Scine {
namespace Sparrow {

namespace nddo {

NDDOElectronicEnergyCalculator::NDDOElectronicEnergyCalculator(const Utils::DensityMatrix& densityMatrix,
                                                               const FockMatrix& fockCalculator,
                                                               const bool& unrestrictedCalculationRunning)
  : densityMatrix_(densityMatrix),
    oneElectronMatrix_(fockCalculator.getOneElectronMatrix()),
    twoElectronMatrix_(fockCalculator.getTwoElectronMatrix()),
    unrestrictedCalculationRunning_(unrestrictedCalculationRunning),
    densityIndependentContributions_(fockCalculator.getDensityIndependentContributions()),
    densityDependentContributions_(fockCalculator.getDensityDependentContributions()) {
}

double NDDOElectronicEnergyCalculator::calculateElectronicEnergy() {
  for (auto const& contribution : densityDependentContributions_) {
    if (contribution->isValid()) {
      contribution->calculate(densityMatrix_, Utils::DerivativeOrder::Zero);
    }
  }
  for (auto const& contribution : densityIndependentContributions_) {
    if (contribution->isValid()) {
      contribution->calculate(densityMatrix_, Utils::DerivativeOrder::Zero);
    }
  }
  if (!unrestrictedCalculationRunning_) {
    return restrictedEnergy();
  }
  return unrestrictedEnergy();
}

double NDDOElectronicEnergyCalculator::restrictedEnergy() {
  auto nAOs = densityMatrix_.restrictedMatrix().rows();
  double electronicEnergy = 0;

  for (unsigned int i = 0; i < nAOs; i++) {
    electronicEnergy += 0.5 * densityMatrix_.restricted(i, i) *
                        (twoElectronMatrix_.getMatrix()(i, i) + 2 * oneElectronMatrix_.getMatrix()(i, i));
    for (unsigned int j = 0; j < i; j++) {
      electronicEnergy += densityMatrix_.restricted(i, j) *
                          (twoElectronMatrix_.getMatrix()(i, j) + 2 * oneElectronMatrix_.getMatrix()(i, j));
    }
  }
  for (auto const& contribution : densityDependentContributions_) {
    electronicEnergy += contribution->getElectronicEnergyContribution();
  }
  for (auto const& contribution : densityIndependentContributions_) {
    electronicEnergy += contribution->getElectronicEnergyContribution();
  }
  return electronicEnergy;
}

double NDDOElectronicEnergyCalculator::unrestrictedEnergy() {
  auto nAOs = densityMatrix_.restrictedMatrix().rows();
  double electronicEnergy = 0;
  for (unsigned int i = 0; i < nAOs; i++) {
    electronicEnergy += 0.5 * (densityMatrix_.restricted(i, i) * (2 * oneElectronMatrix_.getMatrix()(i, i)) +
                               densityMatrix_.alpha(i, i) * twoElectronMatrix_.getAlpha()(i, i) +
                               densityMatrix_.beta(i, i) * twoElectronMatrix_.getBeta()(i, i));
    for (unsigned int j = 0; j < i; j++) {
      electronicEnergy += densityMatrix_.restricted(i, j) * 2 * oneElectronMatrix_.getMatrix()(i, j) +
                          densityMatrix_.alpha(i, j) * twoElectronMatrix_.getAlpha()(i, j) +
                          densityMatrix_.beta(i, j) * twoElectronMatrix_.getBeta()(i, j);
    }
  }
  for (auto const& contribution : densityDependentContributions_) {
    electronicEnergy += contribution->getElectronicEnergyContribution();
  }
  for (auto const& contribution : densityIndependentContributions_) {
    electronicEnergy += contribution->getElectronicEnergyContribution();
  }
  return electronicEnergy;
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
