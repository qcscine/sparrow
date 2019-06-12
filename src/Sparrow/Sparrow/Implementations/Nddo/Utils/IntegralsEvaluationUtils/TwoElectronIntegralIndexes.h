/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_TWOELECTRONINTEGRALINDEXES_H
#define SPARROW_TWOELECTRONINTEGRALINDEXES_H

#include <array>

namespace Scine {
namespace Sparrow {

namespace nddo {

/**
 * @brief This class initializes and stores as a static array the indices of the charge distributions
 *        on one center playing a role in the multipole expansion.
 *
 * The generation of the indices array happens only one time being it declared static.
 * The generated charge distributions are located on one single center, as in the NDDO approximation charge
 * distributions centered on two centers are neglected:
 * \f$ \langle \phi_\mu\phi_\nu | \phi_\\lambda\phi_\sigma \rangle = \vardelta_{IJ}\vardelta_{KL}\langle
 * \chi_\mu^I\chi_\nu^J | \chi_\\lambda^K\chi_\sigma^L\rangle \f$ In total, 40 possible charge distributions are
 * possible: out of a symmetric 9x9 matrix, the diagonal and the terms that give rise to a charge distribution are
 * retained. In the end 40 indices are stored. 5 charge distributions do not give rise to a multipole:
 * - \f$ p_x d_{yz} \f$
 * - \f$ p_y d_{xz} \f$
 * - \f$ p_z d_{xy} \f$
 * - \f$ p_z d_{x^2 - y^2} \f$
 * - \f$ p_xy d_{x^2 - y^2} \f$
 */

class TwoElectronIntegralIndexes {
 public:
  using orb_index_t = int;
  using orbPair_index_t = int;
  /// @brief 9x9 matrix with all the possible orbital combinations
  using Array = std::array<std::array<orbPair_index_t, 9>, 9>;

  /** @brief Convert two orbitals into an orbital pair
   *
   * @param i1 An element of the enum orb_t, describing a single orbital, casted to an integer
   * @param i2 An element of the enum orb_t, describing a single orbital, casted to an integer
   * @return The corresponding element in the array, casted as an integer. The integer 100 corresponds
   *         to an invalid pair (i.e. a pair with no interaction)
   */
  static orbPair_index_t getPairIndex(orb_index_t i1, orb_index_t i2) {
    return getIndex(i1, i2);
  }

  /**
   * @brief Checks if two orbitals form a valid charge distribution
   * @param i1 An element of the enum orb_t, describing a single orbital, casted to an integer
   * @param i2 An element of the enum orb_t, describing a single orbital, casted to an integer
   * @return A boolean indicating if the given pair forms a valid charge distribution
   */
  static bool pairIsInvalid(orb_index_t i1, orb_index_t i2) {
    return (getIndex(i1, i2) == invalidPair);
  }

 private:
  // Creates the index array if it was not yet done and returns the index
  static orbPair_index_t getIndex(orb_index_t i1, orb_index_t i2) {
    static Array pairIndexes = createUniqueIndexes();
    return pairIndexes[i1][i2];
  }

  // Array is created as exclusively invalidPair code except the valid pairs
  // where the casted orbital pair value is given
  static Array createUniqueIndexes();

  // Code for the invalid pair
  static constexpr orbPair_index_t invalidPair = 100;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_TWOELECTRONINTEGRALINDEXES_H