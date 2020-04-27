/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_NDDO_ELECTRONICENERGYCALCULATOR_H
#define SPARROW_NDDO_ELECTRONICENERGYCALCULATOR_H

#include <Utils/Scf/MethodInterfaces/ElectronicEnergyCalculator.h>
#include <memory>
#include <vector>

namespace Scine {

namespace Utils {
class DensityMatrix;
class AdditiveElectronicContribution;
} // namespace Utils

namespace Sparrow {

namespace nddo {
class FockMatrix;
class OneElectronMatrix;
class TwoElectronMatrix;

class NDDOElectronicEnergyCalculator : public Utils::ElectronicEnergyCalculator {
 public:
  NDDOElectronicEnergyCalculator(const Utils::DensityMatrix& densityMatrix, const FockMatrix& fockCalculator,
                                 const bool& unrestrictedCalculationRunning);

  double calculateElectronicEnergy() override;

 private:
  double restrictedEnergy();
  double unrestrictedEnergy();

  const Utils::DensityMatrix& densityMatrix_;
  const OneElectronMatrix& oneElectronMatrix_;
  const TwoElectronMatrix& twoElectronMatrix_;
  const bool& unrestrictedCalculationRunning_;
  const std::vector<std::shared_ptr<Utils::AdditiveElectronicContribution>>&densityIndependentContributions_,
      densityDependentContributions_;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_NDDO_ELECTRONICENERGYCALCULATOR_H
