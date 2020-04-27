/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB_ZEROORDERMATRICESCALCULATOR_H
#define SPARROW_DFTB_ZEROORDERMATRICESCALCULATOR_H

#include "DFTBCommon.h"
#include <Utils/DataStructures/MatrixWithDerivatives.h>
#include <Utils/Math/AutomaticDifferentiation/MethodsTypesHelper.h>
#include <Eigen/Core>

namespace Scine {

namespace Utils {
enum class derivOrder;
class DensityMatrix;
class AtomsOrbitalsIndexes;
} // namespace Utils

namespace Sparrow {

namespace dftb {

/**
 * @brief This class calculates the matrices resulting from the zeroth order expansion of the DFT energy for the DFTB
 *        methods.
 */
class ZeroOrderMatricesCalculator {
 public:
  ZeroOrderMatricesCalculator(const Utils::ElementTypeCollection& elements, const Utils::PositionCollection& positions,
                              const Utils::AtomsOrbitalsIndexes& aoIndexes, const DFTBCommon::AtomicParameterContainer& atomicPar,
                              const DFTBCommon::DiatomicParameterContainer& diatomicPar,
                              const Utils::DensityMatrix& densityMatrix);
  /**
   * @brief Initializes the zeroth order corrected Hamiltonian and the overlap matrices.
   * Furthermore, the one-center blocks are already pre-calculated, as they remain constant.
   */
  void initializeH0S();
  //! @brief Calculate the remaining parts of the Hamiltonian and overlap matrices.
  void constructH0S(Utils::derivOrder order);

  /**
   * @brief Correspond to functions from ElectronicContributionCalculator
   *  Some of these functions do nothing to avoid double initialization (also called for the overlap)
   */
  void initializeFockCalculator() {
  }
  void calculateFockMatrix(Utils::derivOrder /*order*/) {
  }
  // For DFTB0, overlapDerivativeMultiplier is the energy-weighted density matrix.
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::first>& derivatives,
                      const Eigen::MatrixXd& overlapDerivativeMultiplier) const;
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives,
                      const Eigen::MatrixXd& overlapDerivativeMultiplier) const;
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_full>& derivatives,
                      const Eigen::MatrixXd& overlapDerivativeMultiplier) const;

  // Correspond to functions from OverlapCalculator
  void calculateOverlap(Utils::derivOrder highestRequiredOrder);
  const Utils::MatrixWithDerivatives& getOverlap() const;
  const Utils::MatrixWithDerivatives& getZeroOrderHamiltonian() const;
  void resetOverlap();

 private:
  template<Utils::derivOrder O>
  void constructH0S();
  template<Utils::derivOrder O>
  void constructPartOfH0S();
  template<Utils::derivativeType O>
  void addDerivativesImpl(Utils::AutomaticDifferentiation::DerivativeContainerType<O>& derivatives,
                          const Eigen::MatrixXd& overlapDerivativeMultiplier) const;

  Utils::MatrixWithDerivatives zeroOrderHamiltonian_;
  Utils::MatrixWithDerivatives overlap_;
  const Utils::ElementTypeCollection& elements_;
  const Utils::PositionCollection& positions_;
  const Utils::AtomsOrbitalsIndexes& aoIndexes_;
  const DFTBCommon::AtomicParameterContainer& atomicPar_;
  const DFTBCommon::DiatomicParameterContainer& diatomicPar_;
  const Utils::DensityMatrix& densityMatrix_;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB_ZEROORDERMATRICESCALCULATOR_H