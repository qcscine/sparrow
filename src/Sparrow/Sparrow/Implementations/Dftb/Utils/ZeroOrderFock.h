/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB_ZEROORDERFOCK_H
#define SPARROW_DFTB_ZEROORDERFOCK_H

#include <Utils/MethodEssentials/Methods/ElectronicContributionCalculator.h>
#include <Eigen/Core>

namespace Scine {

namespace Utils {
class SingleParticleEnergies;
}

namespace Sparrow {

namespace dftb {
class ZeroOrderMatricesCalculator;

/*!
 * Implementation of FockMatrixCalculator for DFTB0.
 */
class ZeroOrderFock : public Utils::ElectronicContributionCalculator {
 public:
  explicit ZeroOrderFock(ZeroOrderMatricesCalculator& matricesCalculator,
                         const Utils::SingleParticleEnergies& singleParticleEnergies,
                         const Eigen::MatrixXd& energyWeightedDensityMatrix, const int& nElectrons);

  void initialize() override;
  void calculateDensityIndependentPart(Utils::derivOrder order) override;
  void calculateDensityDependentPart(Utils::derivOrder order) override;
  void finalize(Utils::derivOrder /*order*/) override {
  }
  Utils::SpinAdaptedMatrix getMatrix() const override;
  double calculateElectronicEnergy() const override;
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::first>& derivatives) const override;
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives) const override;
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_full>& derivatives) const override;

 private:
  ZeroOrderMatricesCalculator& matricesCalculator_;
  const Utils::SingleParticleEnergies& singleParticleEnergies_;
  const Eigen::MatrixXd& energyWeightedDensityMatrix_;
  const int& nElectrons_;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB_ZEROORDERFOCK_H