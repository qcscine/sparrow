/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_DIPOLEMATRIXCALCULATOR_H
#define SPARROW_DIPOLEMATRIXCALCULATOR_H

#include <Utils/Typenames.h>
#include <exception>

namespace Scine {
namespace Utils {
class DipoleMatrix;
}
namespace Sparrow {

class DipoleMatrixTypeNotAvailableException : public std::exception {
 public:
  const char* what() const noexcept final {
    return "Matrix type not available";
  }
};

/**
 * @brief List of the available calculation methods for calculating the density matrix elements.
 */
enum class IntegralMethod { ObaraSaika, ClosedForm };

/**
 * @brief Interface for the calculation of the dipole matrix in semiempirical methods.
 */
class DipoleMatrixCalculator {
 public:
  virtual ~DipoleMatrixCalculator() = default;
  /**
   * @brief Getter for the dipole matrix in atomic orbitals (AO) basis.
   * @return The dipole matrix in AO basis.
   */
  virtual const Utils::DipoleMatrix& getAODipoleMatrix() const = 0;
  /**
   * @brief Getter for the dipole matrix in moleculat orbital (MO) basis.
   * @return The dipole matrix in MO basis.
   */
  virtual Utils::DipoleMatrix getMODipoleMatrix() const = 0;
  /**
   * @brief Setter for the dipole matrix in AO basis.
   */
  virtual void setAODipoleMatrix(Utils::DipoleMatrix dipoleMatrix) = 0;
  /**
   * @brief Calculates the dipole matrix in AO or MO basis.
   * The method implementation decides if the dipole matrix is calculated in AO or MO basis:
   * while the dipole integral over atomic orbitals is available for NDDO methods, the same does not hold
   * for DFTB methods. For the latter the dipole matrix is calculated directly in MO basis.
   */
  virtual void fillDipoleMatrix(const Eigen::RowVector3d& dipoleEvaluationCoordinate) = 0;

  /**
   * @brief Initialized the dipole matrix calculator and invalidates the current dipole matrix.
   */
  virtual void initialize() = 0;
  /**
   * @brief Setter for the desired dipole matrix integrals calculation method.
   * This is right now only used by the NDDO method.
   */
  virtual void setIntegralMethod(const IntegralMethod& IntegralMethod) = 0;
  /**
   * @brief Getter for the validity status of the dipole matrix.
   */
  virtual bool isValid() const = 0;
  /**
   * @brief Invalidates the underlying dipole matrices and forces a new calculation.
   */
  virtual void invalidate() = 0;
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_DIPOLEMATRIXCALCULATOR_H
