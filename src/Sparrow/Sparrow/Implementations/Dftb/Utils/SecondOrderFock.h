/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB_SECONDORDERFOCK_H
#define SPARROW_DFTB_SECONDORDERFOCK_H

#include "DFTBCommon.h"
#include "ScfFock.h"
#include <Utils/DataStructures/MatrixWithDerivatives.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>

namespace Scine {
namespace Sparrow {

namespace dftb {

/**
 * @brief Implementation of FockMatrixCalculator for DFTB2, the SCC-DFTB.
 *        It calculates the electronic contributions to the energy and their derivatives
 *        with respect to the nuclear cartesian coordinates.
 */
class SecondOrderFock : public ScfFock {
 public:
  //! @brief Constructor calling the ScfFock constructor.
  explicit SecondOrderFock(ZeroOrderMatricesCalculator& matricesCalculator, const Utils::ElementTypeCollection& elements,
                           const Utils::PositionCollection& positions, const DFTBCommon::AtomicParameterContainer& atomicPar,
                           const DFTBCommon::DiatomicParameterContainer& diatomicPar, const Utils::DensityMatrix& densityMatrix,
                           const Eigen::MatrixXd& energyWeightedDensityMatrix, std::vector<double>& atomicCharges,
                           const std::vector<double>& coreCharges, const Utils::AtomsOrbitalsIndexes& aoIndexes,
                           const Eigen::MatrixXd& overlapMatrix, const bool& unrestrictedCalculationRunning);

  //! @brief Initializes the gamma matrix G and the gamma matrix with the derivatives dG.
  void initialize() override;
  //! @brief Sums up the electronic contributions of the zeroth order Hamiltonian, the gamma matrix and the spin
  //! contribution.
  double calculateElectronicEnergy() const override;
  //! @brief adds the derivatives for the first, second atomic and second full types.
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::first>& derivatives) const override;
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives) const override;
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_full>& derivatives) const override;

 protected:
 private:
  /// completes the H matrix by adding the first order correction to H0.
  void completeH() override;
  void constructG(Utils::derivOrder order) override;
  template<Utils::derivOrder O>
  void constructG();
  /*! Return gamma and its derivative(s). */
  template<Utils::derivOrder O>
  Utils::AutomaticDifferentiation::Value1DType<O> gamma(int a, int b) const;
  template<Utils::derivativeType O>
  void addSecondOrderDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<O>& derivatives) const;

  Eigen::MatrixXd G;               // Gamma matrix
  Utils::MatrixWithDerivatives dG; // Derivative of G matrix elements
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB_SECONDORDERFOCK_H