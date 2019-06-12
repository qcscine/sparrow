/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_MULTIPOLEMULTIPOLEINTERACTIONCONTAINER_H
#define SPARROW_MULTIPOLEMULTIPOLEINTERACTIONCONTAINER_H

#include "MultipoleMultipoleInteraction.h"
#include "multipoleTypes.h"
#include <array>

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

/**
 * @brief This class keeps a list of terms of charge-charge-interactions for a pair of multipoles
 *
 * All the possible interactions can be stored in a 13x13 matrix, as there are 13 multipole types:
 * 1 monopole, 3 dipoles, 3 linear quadrupoles, 3 square quadrupoles with charges between the axis, 3 square quadrupoles
 * with charges along the axis.
 */

class MultipoleMultipoleInteractionContainer {
 public:
  //! InteractionLists is the matrix storing the 13x13 possible interactions
  using InteractionLists = std::array<std::array<MultipoleMultipoleInteraction, 13>, 13>;

  /**
   * @brief this function calculates all the possible interactions between the two multipoles m1 and m2
   * @tparam O The order of the derivative required
   * @param m1 The first multipole
   * @param m2 The second multipole
   * @param R The distance between the two multipole centers
   * @param D1 The charge separation in the first multipole
   * @param d2 The charge separation in the second multipole
   * @param squaredRhos The squared Klopman-Ohno parameters
   * @return an Utils::AutomaticDifferentiation::Value1D<O>, i.e. a collection of all the derivative orders up to
   *         the O-th for a 1-dimensional object
   */
  template<Utils::derivOrder O>
  static Utils::AutomaticDifferentiation::Value1DType<O> calculate(multipole_t m1, multipole_t m2, double R, double D1,
                                                                   double d2, double squaredRhos);
  /** @brief The function generated the list of multipole interactions statically, and then returns a
   *         MultipoleMultipoleInteraction object with these two multipoles
   */
  static const MultipoleMultipoleInteraction& get(multipole_t m1, multipole_t m2);

 private:
  // Generates the list of monopoles. Is called to generate a static storage object.
  static InteractionLists generateTermLists();
  // Generates the charge configurations of the multipoles and constructs a MultipoleMultipoleInteraction based on it.
  static MultipoleMultipoleInteraction termList(multipole_t m1, multipole_t m2);
};

inline const MultipoleMultipoleInteraction& MultipoleMultipoleInteractionContainer::get(multipole_t m1, multipole_t m2) {
  static InteractionLists termLists = generateTermLists();
  return termLists[static_cast<int>(m1)][static_cast<int>(m2)];
}

} // namespace multipole

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_MULTIPOLEMULTIPOLEINTERACTIONCONTAINER_H
