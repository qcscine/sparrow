/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB_THIRDORDERFOCK_H
#define SPARROW_DFTB_THIRDORDERFOCK_H

#include "DFTBCommon.h"
#include "ScfFock.h"
#include <Utils/DataStructures/MatrixWithDerivatives.h>
#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Sparrow {

namespace dftb {

/**
 * @brief Implementation of FockMatrixCalculator for DFTB3. It calculates the electronic contributions to the energy
 *        and their derivatives with respect to the nuclear cartesian coordinates.
 */
class ThirdOrderFock : public ScfFock {
 public:
  //! @brief Constructor, calls the ScfFock constructor and sets zeta to 4.0.
  explicit ThirdOrderFock(ZeroOrderMatricesCalculator& matricesCalculator, const Utils::ElementTypeCollection& elements,
                          const Utils::PositionCollection& positions, const DFTBCommon::AtomicParameterContainer& atomicPar,
                          const DFTBCommon::DiatomicParameterContainer& diatomicPar, const Utils::DensityMatrix& densityMatrix,
                          const Eigen::MatrixXd& energyWeightedDensityMatrix, std::vector<double>& atomicCharges,
                          const std::vector<double>& coreCharges, const Utils::AtomsOrbitalsIndexes& aoIndexes,
                          const Eigen::MatrixXd& overlapMatrix, const bool& unrestrictedCalculationRunning);

  //! @brief Calls ScfFock::initialize() and initializes the gamma, Gamma and the matrices containing the derivatives.
  void initialize() override;
  //! @brief Sums up the electronic energy contributions of the zeroth, first and second order Hamiltonian corrections.
  double calculateElectronicEnergy() const override;
  //! @brief calculates automatically the derivatives of the energy with respect to the nuclear cartesian coordinates.
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::first>& derivatives) const override;
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives) const override;
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_full>& derivatives) const override;

 protected:
 private:
  void completeH() override;
  void constructG(Utils::derivOrder order) override;
  template<Utils::derivOrder O>
  void constructG();
  template<Utils::derivOrder O>
  void gammah(int a, int b, Utils::AutomaticDifferentiation::Value1DType<O>& gamma,
              Utils::AutomaticDifferentiation::Value1DType<O>& Gab, Utils::AutomaticDifferentiation::Value1DType<O>& Gba) const;
  template<Utils::derivOrder O>
  void hFactor(double Ua, double Ub, const Utils::AutomaticDifferentiation::Value1DType<O>& R,
               Utils::AutomaticDifferentiation::Value1DType<O>& h, Utils::AutomaticDifferentiation::Value1DType<O>& dhdU) const;
  template<Utils::derivativeType O>
  void addThirdOrderDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<O>& derivatives) const;

  Eigen::MatrixXd g;               // gamma matrix
  Eigen::MatrixXd G;               // Gamma matrix
  Utils::MatrixWithDerivatives dG; // Derivative of G matrix elements
  Utils::MatrixWithDerivatives dg; // Derivative of gamma matrix elements, dgab/dRb
  const double zeta;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB_THIRDORDERFOCK_H