/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DIAGONALPRECONDITIONEREVALUATOR_H
#define SPARROW_DIAGONALPRECONDITIONEREVALUATOR_H
#include "OrderTag.h"
#include <Utils/Math/IterativeDiagonalizer/PreconditionerEvaluator.h>
#include <Utils/Math/IterativeDiagonalizer/SpinAdaptedEigenContainer.h>
namespace Scine {
namespace Utils {
namespace LcaoUtils {
class ElectronicOccupation;
}
class SingleParticleEnergies;
} // namespace Utils
namespace Sparrow {

/**
 * @brief Direct preconditioner calculator.
 * @pre This class expects the occupation to be calculated with the aufbau principle.
 * This class generates the preconditioner vector. For each eigenvalue it creates an element
 * p_k = (H_{k,k} - h_k)^{-1}, where H_{k,k} is the approximated diagonal given by
 * the energy difference between the molecular orbitals involved in the k-th single substitution,
 * and h_k is the k-th guess eigenvalue.
 * If the basis functions are ordered in energetic increasing order, then the preconditioning
 * must also be, and it is ever more effective as more the matrix is diagonally dominant,
 * which is the case normally for linear response problems with an Hamiltonian expanded in
 * basis functions ordered by increasing energy.
 */
class DiagonalPreconditionerEvaluator final : public Utils::PreconditionerEvaluator {
 public:
  /**
   * @brief The constructor generates a preconditioner from an ordered energy difference vector.
   */
  explicit DiagonalPreconditionerEvaluator(const Eigen::VectorXd& energyDifferenceVector);
  /**
   * @brief The constructor generates concatenated energyDifferences from an Unrestricted SpinAdaptedContainer.
   */
  explicit DiagonalPreconditionerEvaluator(
      const Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::VectorXd>& energyDifferenceVector);
  /**
   * @brief The constructor generates energyDifferences from a Restricted SpinAdaptedContainer.
   */
  explicit DiagonalPreconditionerEvaluator(
      const Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd>& energyDifferenceVector);
  /**
   * @brief The constructor orders the energy differences in energetic increasing ordering from an Unrestricted
   * SpinAdaptedContainer.
   * @param OrderTag Tag to give if ordering is needed.
   */
  DiagonalPreconditionerEvaluator(const Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::VectorXd>& energyDifferenceVector,
                                  OrderTag);
  /**
   * @brief The constructor orders the energy differences in energetic increasing ordering from a Restricted
   * SpinAdaptedContainer.
   * @param OrderTag Tag to give if ordering is needed.
   */
  DiagonalPreconditionerEvaluator(const Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd>& energyDifferenceVector,
                                  OrderTag);
  ~DiagonalPreconditionerEvaluator() final = default;
  /**
   * @brief Evaluates the preconditioner vector.
   * The preconditioner elements are p_k = (H_{k,k} - h_k)^{-1},
   * where H_{k,k} is the approximated k-th diagonal element of H,
   * given by the energy difference between the molecular orbitals.
   * H_{k,k} = e_a - e_i
   * @param eigenValues The current guess for the eigenvalues: h_k
   */
  Eigen::VectorXd evaluate(const Eigen::VectorXd& vectorToPrecondition, double eigenvalue) const final;

 private:
  Eigen::VectorXd energyDifferences_;
};
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_UNRESTRICTEDDIAGONALPRECONDITIONEREVALUATOR_H
