/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "TwoElectronIntegralIndexes.h"
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/GeneralTypes.h>

namespace Scine {
namespace Sparrow {

namespace nddo {
using namespace GeneralTypes;

TwoElectronIntegralIndexes::Array TwoElectronIntegralIndexes::createUniqueIndexes() {
  Array pairIndexes;

  // Initialize pairIndexes to invalidPair
  for (auto& row : pairIndexes) {
    for (auto& col : row) {
      col = invalidPair;
    }
  }

  pairIndexes[static_cast<orb_index_t>(orb_t::s)][static_cast<orb_index_t>(orb_t::s)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::s_s);
  pairIndexes[static_cast<orb_index_t>(orb_t::s)][static_cast<orb_index_t>(orb_t::z)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::s_z);
  pairIndexes[static_cast<orb_index_t>(orb_t::s)][static_cast<orb_index_t>(orb_t::x)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::s_x);
  pairIndexes[static_cast<orb_index_t>(orb_t::s)][static_cast<orb_index_t>(orb_t::y)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::s_y);
  pairIndexes[static_cast<orb_index_t>(orb_t::z)][static_cast<orb_index_t>(orb_t::z)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::z_z);
  pairIndexes[static_cast<orb_index_t>(orb_t::x)][static_cast<orb_index_t>(orb_t::z)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::x_z);
  pairIndexes[static_cast<orb_index_t>(orb_t::y)][static_cast<orb_index_t>(orb_t::z)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::y_z);
  pairIndexes[static_cast<orb_index_t>(orb_t::x)][static_cast<orb_index_t>(orb_t::x)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::x_x);
  pairIndexes[static_cast<orb_index_t>(orb_t::x)][static_cast<orb_index_t>(orb_t::y)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::x_y);
  pairIndexes[static_cast<orb_index_t>(orb_t::y)][static_cast<orb_index_t>(orb_t::y)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::y_y);
  pairIndexes[static_cast<orb_index_t>(orb_t::s)][static_cast<orb_index_t>(orb_t::z2)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::s_z2);
  pairIndexes[static_cast<orb_index_t>(orb_t::s)][static_cast<orb_index_t>(orb_t::xz)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::s_xz);
  pairIndexes[static_cast<orb_index_t>(orb_t::s)][static_cast<orb_index_t>(orb_t::yz)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::s_yz);
  pairIndexes[static_cast<orb_index_t>(orb_t::s)][static_cast<orb_index_t>(orb_t::x2y2)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::s_x2y2);
  pairIndexes[static_cast<orb_index_t>(orb_t::s)][static_cast<orb_index_t>(orb_t::xy)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::s_xy);
  pairIndexes[static_cast<orb_index_t>(orb_t::z)][static_cast<orb_index_t>(orb_t::z2)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::z_z2);
  pairIndexes[static_cast<orb_index_t>(orb_t::z)][static_cast<orb_index_t>(orb_t::xz)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::z_xz);
  pairIndexes[static_cast<orb_index_t>(orb_t::z)][static_cast<orb_index_t>(orb_t::yz)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::z_yz);
  pairIndexes[static_cast<orb_index_t>(orb_t::x)][static_cast<orb_index_t>(orb_t::z2)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::x_z2);
  pairIndexes[static_cast<orb_index_t>(orb_t::x)][static_cast<orb_index_t>(orb_t::xz)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::x_xz);
  pairIndexes[static_cast<orb_index_t>(orb_t::x)][static_cast<orb_index_t>(orb_t::x2y2)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::x_x2y2);
  pairIndexes[static_cast<orb_index_t>(orb_t::x)][static_cast<orb_index_t>(orb_t::xy)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::x_xy);
  pairIndexes[static_cast<orb_index_t>(orb_t::y)][static_cast<orb_index_t>(orb_t::z2)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::y_z2);
  pairIndexes[static_cast<orb_index_t>(orb_t::y)][static_cast<orb_index_t>(orb_t::yz)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::y_yz);
  pairIndexes[static_cast<orb_index_t>(orb_t::y)][static_cast<orb_index_t>(orb_t::x2y2)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::y_x2y2);
  pairIndexes[static_cast<orb_index_t>(orb_t::y)][static_cast<orb_index_t>(orb_t::xy)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::y_xy);
  pairIndexes[static_cast<orb_index_t>(orb_t::z2)][static_cast<orb_index_t>(orb_t::z2)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::z2_z2);
  pairIndexes[static_cast<orb_index_t>(orb_t::z2)][static_cast<orb_index_t>(orb_t::xz)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::z2_xz);
  pairIndexes[static_cast<orb_index_t>(orb_t::z2)][static_cast<orb_index_t>(orb_t::yz)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::z2_yz);
  pairIndexes[static_cast<orb_index_t>(orb_t::z2)][static_cast<orb_index_t>(orb_t::x2y2)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::z2_x2y2);
  pairIndexes[static_cast<orb_index_t>(orb_t::z2)][static_cast<orb_index_t>(orb_t::xy)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::z2_xy);
  pairIndexes[static_cast<orb_index_t>(orb_t::xz)][static_cast<orb_index_t>(orb_t::xz)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::xz_xz);
  pairIndexes[static_cast<orb_index_t>(orb_t::xz)][static_cast<orb_index_t>(orb_t::yz)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::xz_yz);
  pairIndexes[static_cast<orb_index_t>(orb_t::xz)][static_cast<orb_index_t>(orb_t::x2y2)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::xz_x2y2);
  pairIndexes[static_cast<orb_index_t>(orb_t::xz)][static_cast<orb_index_t>(orb_t::xy)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::xz_xy);
  pairIndexes[static_cast<orb_index_t>(orb_t::yz)][static_cast<orb_index_t>(orb_t::yz)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::yz_yz);
  pairIndexes[static_cast<orb_index_t>(orb_t::yz)][static_cast<orb_index_t>(orb_t::x2y2)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::yz_x2y2);
  pairIndexes[static_cast<orb_index_t>(orb_t::yz)][static_cast<orb_index_t>(orb_t::xy)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::yz_xy);
  pairIndexes[static_cast<orb_index_t>(orb_t::x2y2)][static_cast<orb_index_t>(orb_t::x2y2)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::x2y2_x2y2);
  pairIndexes[static_cast<orb_index_t>(orb_t::xy)][static_cast<orb_index_t>(orb_t::xy)] =
      static_cast<orbPair_index_t>(twoElIntegral_t::xy_xy);

  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      if (pairIndexes[i][j] != invalidPair) {
        pairIndexes[j][i] = pairIndexes[i][j];
      }
    }
  }

  return pairIndexes;
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
