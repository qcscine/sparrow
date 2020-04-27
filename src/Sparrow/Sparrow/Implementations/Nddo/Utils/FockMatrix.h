/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_NDDO_FOCKMATRIX_H
#define SPARROW_NDDO_FOCKMATRIX_H

#include "OneElectronMatrix.h"
#include "TwoElectronMatrix.h"
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/TwoCenterIntegralContainer.h>
#include <Utils/Math/AutomaticDifferentiation/MethodsTypesHelper.h>
#include <Utils/Scf/MethodInterfaces/ElectronicContributionCalculator.h>
#include <memory>

namespace Scine {

namespace Utils {
class AtomsOrbitalsIndexes;
class DensityMatrix;
class OverlapCalculator;
class ElectronicEnergyCalculator;
} // namespace Utils

namespace Sparrow {

namespace nddo {

class FockMatrix : public Utils::ElectronicContributionCalculator {
 public:
  FockMatrix(const Utils::ElementTypeCollection& elements, const Utils::PositionCollection& positions,
             const Utils::DensityMatrix& densityMatrix, const OneCenterIntegralContainer& oneCIntegrals,
             const ElementParameters& elementPar, const Utils::AtomsOrbitalsIndexes& aoIndexes,
             const Utils::OverlapCalculator& overlapCalculator, const bool& unrestrictedCalculationRunning);

  void initialize() override;
  void calculateDensityIndependentPart(Utils::derivOrder order) override;
  void calculateDensityDependentPart(Utils::derivOrder order) override;
  void finalize(Utils::derivOrder order) override;
  Utils::SpinAdaptedMatrix getMatrix() const override;
  double calculateElectronicEnergy() const override;
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::first>& derivatives) const override;
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives) const override;
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_full>& derivatives) const override;

  const OneElectronMatrix& getOneElectronMatrix() const;
  const TwoElectronMatrix& getTwoElectronMatrix() const;
  const std::vector<std::shared_ptr<Utils::AdditiveElectronicContribution>>& getDensityDependentContributions() const;
  const std::vector<std::shared_ptr<Utils::AdditiveElectronicContribution>>& getDensityIndependentContributions() const;

  /**
   * This function adds an additive electronic contribution to the
   * Hamiltonian that will be evaluated each SCF iteration.
   */
  void addDensityDependentElectronicContribution(std::shared_ptr<Utils::AdditiveElectronicContribution> contribution) final;
  /**
   * This function adds an additive electronic contribution to the Hamiltonian
   * that will be evaluated once per single-point calculation.
   */
  void addDensityIndependentElectronicContribution(std::shared_ptr<Utils::AdditiveElectronicContribution> contribution) final;
  void clearElectronicContributions();
  void eraseElectronicContribution(std::shared_ptr<Utils::AdditiveElectronicContribution> contribution);

 private:
  template<Utils::derivativeType O>
  void addDerivativesImpl(Utils::AutomaticDifferentiation::DerivativeContainerType<O>& derivatives) const;

  TwoCenterIntegralContainer twoCenterIntegrals_;
  OneElectronMatrix F1_;
  TwoElectronMatrix F2_;
  const Utils::OverlapCalculator& overlapCalculator_;
  const bool& unrestrictedCalculationRunning_;
  std::unique_ptr<Utils::ElectronicEnergyCalculator> electronicEnergyCalculator_;
  std::vector<std::shared_ptr<Utils::AdditiveElectronicContribution>> densityDependentContributions_,
      densityIndependentContributions_;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_NDDO_FOCKMATRIX_H
