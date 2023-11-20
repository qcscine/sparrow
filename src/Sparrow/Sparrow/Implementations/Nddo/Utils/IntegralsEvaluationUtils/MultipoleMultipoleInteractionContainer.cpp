/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MultipoleMultipoleInteractionContainer.h"
#include "ChargesInMultipoles.h"
#include "MMTermCreator.h"

namespace Scine {
namespace Sparrow {

using namespace Utils::AutomaticDifferentiation;

namespace nddo {

namespace multipole {

MultipoleMultipoleInteractionContainer::InteractionLists MultipoleMultipoleInteractionContainer::generateTermLists() {
  InteractionLists termLists;
  // Created the 13x13 matrix corresponding to the possible multipole-multipole interactions
  for (int i = 0; i < 13; i++) {
    for (int j = 0; j < 13; j++) { // TODO: symmetry ?
      termLists[i][j] = termList(static_cast<Multipole>(i), static_cast<Multipole>(j));
    }
  }

  return termLists;
}

MultipoleMultipoleInteraction MultipoleMultipoleInteractionContainer::termList(Multipole m1, Multipole m2) {
  MMTermCreator creator;
  const auto& charges1 = ChargesInMultipoles::getChargeConfiguration(m1);
  const auto& charges2 = ChargesInMultipoles::getChargeConfiguration(m2);

  for (const auto& c1 : charges1) {
    for (const auto& c2 : charges2) {
      creator.add(MultipoleChargePair(c1, c2));
    }
  }

  return creator.computeList();
}

template<Utils::DerivativeOrder O>
Value1DType<O> MultipoleMultipoleInteractionContainer::calculate(Multipole m1, Multipole m2, double R, double D1,
                                                                 double d2, double squaredRhos) {
  return get(m1, m2).calculate<O>(R, D1, d2, squaredRhos);
}

template Value1DType<Utils::DerivativeOrder::Zero>
MultipoleMultipoleInteractionContainer::calculate<Utils::DerivativeOrder::Zero>(Multipole m1, Multipole m2, double R,
                                                                                double D1, double d2, double squaredRhos);
template Value1DType<Utils::DerivativeOrder::One>
MultipoleMultipoleInteractionContainer::calculate<Utils::DerivativeOrder::One>(Multipole m1, Multipole m2, double R,
                                                                               double D1, double d2, double squaredRhos);
template Value1DType<Utils::DerivativeOrder::Two>
MultipoleMultipoleInteractionContainer::calculate<Utils::DerivativeOrder::Two>(Multipole m1, Multipole m2, double R,
                                                                               double D1, double d2, double squaredRhos);
} // namespace multipole

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
