/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_ATOMPAIRDIPOLE_H
#define SPARROW_ATOMPAIRDIPOLE_H

#include <Utils/Typenames.h>
#include <Eigen/Core>
#include <array>

namespace Scine {

namespace Utils {
class AtomicGtos;
class DipoleMatrix;
} // namespace Utils

namespace Sparrow {
enum class IntegralMethod;

/**
 * @brief Class responsible for calculating a block of the dipole matrix
 */
class AtomPairDipole {
 public:
  /**
   * @brief Calculates a block of the dipole matrix between two atoms.
   * @param dipoleMatrix The dipole matrix given as reference.
   * @param startOfAtomA The index corresponding to the first atomic orbital of the atom A.
   * @param startOfAtomB The index corresponding to the first atomic orbital of the atom B.
   * @param method Decides which method to use for the calculation of the integrals.
   * @param gtosA The GTO expansion on the atom A.
   * @param gtosB The GTO expansion on the atom B.
   * @param Ra Positions of the nucleus A.
   * @param Rb Positions of the nucleus B.
   * @param Rab Vectorial separation between A and B.
   * @param dipoleEvaluationCoordinate Decides where the dipole has to be calculated from.
   */
  static void fillAtomPairDipoleBlock(Utils::DipoleMatrix& dipoleMatrix, int startOfAtomA, int startOfAtomB,
                                      const IntegralMethod& method, const Utils::AtomicGtos& gtosA,
                                      const Utils::AtomicGtos& gtosB, const Eigen::RowVector3d& Ra,
                                      const Eigen::RowVector3d& Rb, const Eigen::RowVector3d& Rab,
                                      const Eigen::RowVector3d& dipoleEvaluationCoordinate);
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_ATOMPAIRDIPOLE_H
