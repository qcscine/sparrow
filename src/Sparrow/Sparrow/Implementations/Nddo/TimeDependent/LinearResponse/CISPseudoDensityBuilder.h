/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_CISPSEUDODENSITYBUILDER_H
#define SPARROW_CISPSEUDODENSITYBUILDER_H

#include "CISData.h"
#include <Sparrow/Implementations/TimeDependent/TimeDependentUtils.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
namespace Scine {
namespace Sparrow {

/**

 * @brief This class evaluates pseudo-density matrices for the construction of the sigma-vectors in the Davidson-Liu
 * algorithm. The calculation is performed according to J. B. Foresman, M. Head-Gordon, J. A. Pople and M. J. Frisch; J.
 * Phys. Chem. 1992, 96, 1, 135-149.
 *
 * \f$P\f$: pseudo-density matrix
 *
 * \f$C\f$: coefficient matrix with either only occupied or virtual orbitals
 *
 * \f$c\f$: guess vector
 *
 * single element of \f$P\f$: \f$ P_{\lambda \sigma} = \sum_{j b}{C_{\lambda j}}^* c_{j b}  C_{\sigma b} \f$
 *
 * Whole \f$P\f$: \f$  P = C^T c~ C \f$
 *
 * The coefficient matrices are built upon construction of the class, since the full coefficient matrix does not change
 * as well as the occupation.
 *
 * @tparam restrictedness of the reference calculation as template argument (either restricted or unrestricted)
 * @class CISPseudoDensityBuilder
 * @file CISPseudoDensityBuilder.h
 */
template<Utils::Reference restrictedness>
class CISPseudoDensityBuilder {
 public:
  using SAMType = Utils::SpinAdaptedContainer<restrictedness, Eigen::MatrixXd>;
  /**
   * @brief Constructs the class by building the virtual and occupied blocks from the full coefficient matrix. This is
   * done by initializing two matrices to be filled by construct orbitals together with the full coefficient matrix and
   * a vector with the indexes of the occupied orbitals. The constructor is template-specified to call
   * construct-Orbitals once for the restricted and twice for the unrestricted reference (alpha and beta).
   * @param molecularOrbitals
   * @param occupation
   */
  CISPseudoDensityBuilder(const Utils::MolecularOrbitals& molecularOrbitals,
                          const Utils::LcaoUtils::ElectronicOccupation& occupation);
  /**
   * @brief Destructor
   */
  ~CISPseudoDensityBuilder() = default;
  /**
   * @brief Constructs the virtual and occupied blocks from the full coefficient matrix.
   * The blocks are constructed from the full coefficient matrix (allMOs) and the index vector (filledOrbitals), as well
   * as two already initialized matrices for the virtual and occupied coefficients.
   * @param occupiedOrbitals
   * @param virtualOrbitals
   * @param allMOs
   * @param filledOrbitals
   */
  void constructOrbitals(Eigen::MatrixXd& occupiedOrbitals, Eigen::MatrixXd& virtualOrbitals,
                         const Eigen::MatrixXd& allMOs, const std::vector<int>& filledOrbitals);
  /**
   *@brief Evaluates and returns the pseudo-density matrix
   * Calculates the pseudo-density matrix with one guessVector from the Davidson-Liu algorithm. This function calls
   * \ref mapAndMultiply(const Eigen::VectorXd &, const Eigen::MatrixXd &, const Eigen::MatrixXd &) const mapAndMultiply
   *for the calculation of the restricted, alpha and beta variants.
   *
   * @param guessVector
   * @return pseudoDensityMatrix
   */
  const SAMType getPseudoDensityMatrix(const Utils::SpinAdaptedContainer<restrictedness, Eigen::VectorXd>& guessVector) const;
  /**
   * @brief Getter for occupied block of coefficient matrix
   * @return occupiedOrbitals
   */
  SAMType getOccupiedOrbitals() const;
  /**
   * @brief Getter for virtual block of coefficient matrix
   * @return virtualOrbitals
   */
  SAMType getVirtualOrbitals() const;
  /**
   * @brief Perform a matrix multiplication of the coefficient blocks with a mapped guessVector.
   * @param Vector
   * @param virtualMOs
   * @param occupiedMOs
   * @return Eigen::MatrixXd
   */
  Eigen::MatrixXd mapAndMultiply(const Eigen::VectorXd& Vector, const Eigen::MatrixXd& virtualMOs,
                                 const Eigen::MatrixXd& occupiedMOs) const;

 private:
  SAMType occupiedMolecularOrbitals_;
  SAMType virtualMolecularOrbitals_;
};

template<Utils::Reference restrictedness>
Utils::SpinAdaptedContainer<restrictedness, Eigen::MatrixXd> CISPseudoDensityBuilder<restrictedness>::getOccupiedOrbitals() const {
  return occupiedMolecularOrbitals_;
}

template<Utils::Reference restrictedness>
typename CISPseudoDensityBuilder<restrictedness>::SAMType CISPseudoDensityBuilder<restrictedness>::getVirtualOrbitals() const {
  return virtualMolecularOrbitals_;
}

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_CISPSEUDODENSITYBUILDER_H
