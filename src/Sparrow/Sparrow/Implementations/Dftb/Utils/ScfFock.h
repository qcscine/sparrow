/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB_SCFFOCK_H
#define SPARROW_DFTB_SCFFOCK_H

#include "DFTBCommon.h"
#include "SDFTB.h"
#include <Utils/MethodEssentials/Methods/ElectronicContributionCalculator.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Utils {
class AtomsOrbitalsIndexes;
class DensityMatrix;
} // namespace Utils

namespace Sparrow {

namespace dftb {
class ZeroOrderMatricesCalculator;

/*!
 * Implementation of FockMatrixCalculator for SCF-type DFTB methods.
 */
class ScfFock : public Utils::ElectronicContributionCalculator {
 public:
  explicit ScfFock(ZeroOrderMatricesCalculator& matricesCalculator, const Utils::ElementTypeCollection& elements,
                   const Utils::PositionCollection& positions, const DFTBCommon::AtomicParameterContainer& atomicPar,
                   const DFTBCommon::DiatomicParameterContainer& diatomicPar, const Utils::DensityMatrix& densityMatrix,
                   const Eigen::MatrixXd& energyWeightedDensityMatrix, std::vector<double>& atomicCharges,
                   const std::vector<double>& coreCharges, const Utils::AtomsOrbitalsIndexes& aoIndexes,
                   const Eigen::MatrixXd& overlapMatrix, const bool& unrestrictedCalculationRunning);

  void initialize() override;
  void calculateDensityDependentPart(Utils::derivOrder order) override;
  void calculateDensityIndependentPart(Utils::derivOrder order) override;
  Utils::SpinAdaptedMatrix getMatrix() const override;
  void finalize(Utils::derivOrder order) override;

 protected:
  int getNumberAtoms() const;
  void populationAnalysis();

  ZeroOrderMatricesCalculator& zeroOrderMatricesCalculator_;
  const Utils::ElementTypeCollection& elements_;
  const Utils::PositionCollection& positions_;
  const DFTBCommon::AtomicParameterContainer& atomicPar_;
  const DFTBCommon::DiatomicParameterContainer& diatomicPar_;
  const Utils::DensityMatrix& densityMatrix_;
  const Eigen::MatrixXd& energyWeightedDensityMatrix_;
  std::vector<double>& atomicCharges_;
  const std::vector<double>& coreCharges_;
  const Utils::AtomsOrbitalsIndexes& aoIndexes_;
  const Eigen::MatrixXd& overlapMatrix_;
  const bool& unrestrictedCalculationRunning_;

  SDFTB spinDFTB;
  Eigen::MatrixXd HXoverS_;
  Eigen::MatrixXd H0_;
  Eigen::MatrixXd correctionToFock;

 private:
  virtual void completeH() = 0;
  virtual void constructG(Utils::derivOrder order) = 0;
};

inline int ScfFock::getNumberAtoms() const {
  return elements_.size();
}

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB_SCFFOCK_H