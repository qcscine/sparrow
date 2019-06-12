/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Global2c2eTerms.h"
#include "zeroLocal2c2eIntegrals.h"
#include <array>

namespace Scine {
namespace Sparrow {

namespace nddo {
using namespace GeneralTypes;

namespace multipole {

Global2c2eTerms::RotationTermsArray Global2c2eTerms::createRotationTerms() {
  RotationTermsArray expressions(40, std::vector<RotationTerms>(40));
  // no s orbital (so not 9 orbitals, just 8)
  std::array<int, 8> i{};
  for (i[0] = 0; i[0] < 9; i[0]++) {
    for (i[1] = i[0]; i[1] < 9; i[1]++) {
      for (i[2] = 0; i[2] < 9; i[2]++) {
        for (i[3] = i[2]; i[3] < 9; i[3]++) {
          for (i[4] = 0; i[4] < 9; i[4]++) {
            for (i[5] = 0; i[5] < 9; i[5]++) {
              for (i[6] = 0; i[6] < 9; i[6]++) {
                for (i[7] = 0; i[7] < 9; i[7]++) {
                  createTerm(expressions, i);
                }
              }
            }
          }
        }
      }
    }
  }

  return expressions;
}

void Global2c2eTerms::createTerm(RotationTermsArray& expressions, std::array<int, 8> i) {
  std::array<int, 8> l{};
  std::array<orbPair_index_t, 4> pair{};
  std::array<orb_t, 8> o{};

  for (int ind = 0; ind < 4; ind++) {
    if (!compatibleOrbitals(i[ind], i[ind + 4]))
      return;

    pair[ind] = TwoElectronIntegralIndexes::getPairIndex(i[2 * ind], i[2 * ind + 1]);
    if (TwoElectronIntegralIndexes::pairIsInvalid(i[2 * ind], i[2 * ind + 1]))
      return;
  }

  if (ZeroLocal2c2eIntegrals::isZero(static_cast<twoElIntegral_t>(pair[2]), static_cast<twoElIntegral_t>(pair[3])))
    return;

  for (int ind = 0; ind < 8; ind++) {
    o[ind] = static_cast<orb_t>(i[ind]);
    l[ind] = orbitalQN(o[ind]);
  }

  // Rotations with factor 0
  for (int ind = 0; ind < 4; ind++) {
    if (l[ind] == 1 && o[ind] == z && o[ind + 4] == y)
      return;
    if (l[ind] == 2 && o[ind] == z2 && o[ind + 4] == yz)
      return;
    if (l[ind] == 2 && o[ind] == z2 && o[ind + 4] == xy)
      return;
  }

  expressions[pair[0]][pair[1]].emplace_back(pair[2], pair[3], getRotPairType(o[0], o[4]), getRotPairType(o[1], o[5]),
                                             getRotPairType(o[2], o[6]), getRotPairType(o[3], o[7]));
}

bool Global2c2eTerms::compatibleOrbitals(int a, int b) {
  if (a == 0 && b == 0)
    return true;
  if (a >= 1 && a <= 3 && b >= 1 && b <= 3)
    return true;
  if (a > 3 && b > 3)
    return true;

  return false;
}

} // namespace multipole

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
