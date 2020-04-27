/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB_ZEROORDERFOCK_H
#define SPARROW_DFTB_ZEROORDERFOCK_H

#include <Utils/Scf/MethodInterfaces/ElectronicContributionCalculator.h>
#include <Eigen/Core>

namespace Scine {

namespace Utils {
class SingleParticleEnergies;
class AdditiveElectronicContribution;
} // namespace Utils

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
  /**
   * This function adds an additive electronic contribution to the
   * Hamiltonian that will be evaluated each SCF iteration.
   * At zero order no SCF is done. This just calls the density independent version
   */
  void addDensityDependentElectronicContribution(std::shared_ptr<Utils::AdditiveElectronicContribution> contribution) final;
  /**
   * This function adds an additive electronic contribution to the Hamiltonian
   * that will be evaluated once per single-point calculation.
   */
  void addDensityIndependentElectronicContribution(std::shared_ptr<Utils::AdditiveElectronicContribution> contribution) final;

 private:
  ZeroOrderMatricesCalculator& matricesCalculator_;
  const Utils::SingleParticleEnergies& singleParticleEnergies_;
  const Eigen::MatrixXd& energyWeightedDensityMatrix_;
  const int& nElectrons_;
  std::vector<std::shared_ptr<Utils::AdditiveElectronicContribution>> densityIndependentContributions_;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB_ZEROORDERFOCK_H