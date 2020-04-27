/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTBDIPOLEMATRIXCALCULATOR_H
#define SPARROW_DFTBDIPOLEMATRIXCALCULATOR_H

#include <Sparrow/Implementations/DipoleMatrixCalculator.h>
#include <Utils/DataStructures/DipoleMatrix.h>
#include <Utils/Typenames.h>
#include <memory>

namespace Scine {

namespace Utils {
class MolecularOrbitals;
class AtomsOrbitalsIndexes;
} // namespace Utils

namespace Sparrow {

/**
 * @brief Class responsible for the calculation of the dipole matrix in MO basis for DFTB methods.
 * This class calculates the dipole matrix in MO basis according to
 * R. Rüger, E. van Lenthe, T. Heine and L. Visscher,
 * Tight-Binding Approximations to Time-Dependent Density Functional Theory
 * - a fast approach for the calculation of electronically excited states,
 * J. Chem. Phys. 144, 184103, 2016, https://doi.org/10.1063/1.4948647
 *
 * It works with the DFTBO0, DFTB2 and DFTB3 methods in restricted and
 * unrestricted formalism.
 */
template<class DFTBMethod>
class DFTBDipoleMatrixCalculator : public DipoleMatrixCalculator {
 public:
  /*
   * @brief Factory method for the DFTBDipoleMatrixCalculator.
   * This returns a unique pointer to a DFTBDipoleMatrixCalculator class.
   */
  static std::unique_ptr<DFTBDipoleMatrixCalculator<DFTBMethod>> create(DFTBMethod& method);
  //! @brief Virtual destructor.
  ~DFTBDipoleMatrixCalculator() final;
  //! @brief Getter for the AO dipole matrix. Throws an exception.
  const Utils::DipoleMatrix& getAODipoleMatrix() const final;
  //! @brief Getter for the MO dipole matrix. Returns an rvalue.
  Utils::DipoleMatrix getMODipoleMatrix() const final;
  //! @brief Setter for the AO dipole matrix. Here it throws an exception.
  void setAODipoleMatrix(Utils::DipoleMatrix dipoleMatrix) final;
  //! @brief Calculates the MO dipole matrix for the DFTB methods.
  void fillDipoleMatrix(const Eigen::RowVector3d& dipoleEvaluationCoordinate) final;

  //! @brief Initialize the underlying dipole matrix.
  void initialize() final;
  //! @brief This does nothing in DFTB.
  void setIntegralMethod(const IntegralMethod& IntegralMethod) final;
  //! @brief Invalidates the underlying dipole matrices and forces a new calculation.
  void invalidate() final;
  //! @brief Checks the validity of the underlying dipole matrix.
  bool isValid() const final;

 private:
  explicit DFTBDipoleMatrixCalculator(DFTBMethod& method);
  void calculateRestrictedTransitionChargeMatrices(const Utils::PositionCollection& positions);
  void calculateUnrestrictedTransitionChargeMatrices(const Utils::PositionCollection& positions);
  bool valid_{false};
  Utils::DipoleMatrix dipoleMatrixMO_;
  const DFTBMethod& method_;

  const Utils::PositionCollection positions_;
  const Utils::MolecularOrbitals& coefficientMatrix_;
  const Eigen::MatrixXd& overlapMatrix_;
  const Utils::AtomsOrbitalsIndexes& aoIndex_;

  std::vector<Eigen::MatrixXd> transitionChargeMatrices_;
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_DFTBDIPOLEMATRIXCALCULATOR_H
